# make_outputs.py
# Put this file in your project root (same level as Taylor.py, error.py, sweep.py, etc.)
# Update the three CSV filenames below if yours are different.

import os
import math
import time
import gc
import statistics
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from pathlib import Path
from scipy.stats import norm
from Taylor import atm_coeff  # must exist in your project

from error import errors_at_x

def build_baseline_sweep_df(S, r, T, sigma, x_min=-0.05, x_max=0.05, step=0.001, save_path=None):
    # get coefficients (atm_coeff might return 4 or 5 values)
    coeff = atm_coeff(S, r, T, sigma)
    if len(coeff) == 4:
        C0, c1, c2, c3 = coeff
    else:
        C0, c1, c2, c3, _ = coeff

    xs = np.arange(x_min, x_max + step, step, dtype=np.float64)

    rows = []
    for x in xs:
        rows.append(errors_at_x(S, r, T, sigma, float(x), C0, c1, c2, c3))

    df = pd.DataFrame(rows)

    # normalize gamma column name if your errors_at_x uses gamma_error but it's actually relative
    if "gamma_rel_error" not in df.columns and "gamma_error" in df.columns:
        df = df.rename(columns={"gamma_error": "gamma_rel_error"})

    if save_path is not None:
        df.to_csv(save_path, index=False)

    return df

# -----------------------------
# Utility: timing
# -----------------------------
def time_one(fn, *args, repeats=5, warmup=1, mode="best"):
    gc_was_enabled = gc.isenabled()
    gc.disable()
    try:
        for _ in range(warmup):
            fn(*args)

        times = []
        for _ in range(repeats):
            t0 = time.perf_counter()
            fn(*args)
            times.append(time.perf_counter() - t0)

        return float(statistics.median(times)) if mode == "median" else float(min(times))
    finally:
        if gc_was_enabled:
            gc.enable()


# -----------------------------
# Vectorized BS and T3 (for speed benchmarks)
# -----------------------------
def bs_vec(S, Ks, r, T, sigma):
    a = sigma * np.sqrt(T)
    d1 = (np.log(S / Ks) + (r + 0.5 * sigma**2) * T) / a
    d2 = d1 - a
    C = S * norm.cdf(d1) - Ks * np.exp(-r * T) * norm.cdf(d2)
    Delta = norm.cdf(d1)
    Gamma = norm.pdf(d1) / (S * a)
    return C, Delta, Gamma

def t3_vec(xs, S, C0, c1, c2, c3):
    C = C0 + xs * (c1 + xs * (0.5 * c2 + xs * (c3 / 6.0)))
    Delta = (c1 + c2 * xs + 0.5 * c3 * xs * xs) / S
    Gamma = ((c2 - c1) + (c3 - c2) * xs - 0.5 * c3 * xs * xs) / (S * S)
    return C, Delta, Gamma


# -----------------------------
# Plot and summary functions (your originals)
# -----------------------------
def plot_error_vs_x(df, x_star, ycol, ybound, title, ylabel, outpath):
    d = df.sort_values("Log-moneyness")
    x = d["Log-moneyness"].to_numpy()
    y = d[ycol].to_numpy()

    plt.figure()
    plt.plot(x, y)
    plt.axhline(ybound)
    plt.axvline(-x_star)
    plt.axvline(+x_star)
    plt.title(title)
    plt.xlabel("log-moneyness x")
    plt.ylabel(ylabel)
    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    plt.close()

def heatmap_xstar(grid_df, outpath):
    pivot = grid_df.pivot(index="T", columns="sigma", values="x_star")
    plt.figure()
    plt.imshow(pivot.to_numpy(), aspect="auto", origin="lower")
    plt.title("Validated band radius x* across (T, sigma)")
    plt.xlabel("sigma")
    plt.ylabel("T")
    plt.xticks(range(len(pivot.columns)), [str(s) for s in pivot.columns])
    plt.yticks(range(len(pivot.index)), [str(t) for t in pivot.index])
    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    plt.close()

def plot_speed(speed_df, outpath_times, outpath_speedup):
    d = speed_df.sort_values("N")

    plt.figure()
    plt.plot(d["N"], d["t_bs_sec"])
    plt.plot(d["N"], d["t_t3_sec"])
    plt.plot(d["N"], d["speedup"], marker="o")
    plt.title("Wall-clock time vs batch size")
    plt.xlabel("N")
    plt.ylabel("seconds")
    plt.yscale("log")
    plt.tight_layout()
    plt.savefig(outpath_times, dpi=200)
    plt.close()

    plt.figure()
    plt.plot(d["N"], d["speedup"])
    plt.title("Speedup (BS time / Taylor time)")
    plt.xlabel("N")
    plt.ylabel("speedup")
    plt.tight_layout()
    plt.savefig(outpath_speedup, dpi=200)
    plt.close()

def write_summary_txt(path, params, thresholds, band_stats, speed_stats):
    p = Path(path)
    lines = []
    lines.append("RESULTS SUMMARY\n")
    lines.append("Parameters:\n")
    for k, v in params.items():
        lines.append(f"  {k}: {v}\n")

    lines.append("\nThresholds:\n")
    for k, v in thresholds.items():
        lines.append(f"  {k}: {v}\n")

    lines.append("\nValidated band (baseline):\n")
    for k, v in band_stats.items():
        lines.append(f"  {k}: {v}\n")

    lines.append("\nSpeed benchmark:\n")
    for row in speed_stats:
        lines.append(
            f"  N={row['N']:,}: BS={row['t_bs_sec']:.4f}s, "
            f"T3={row['t_t3_sec']:.4f}s, speedup={row['speedup']:.2f}x\n"
        )

    p.write_text("".join(lines), encoding="utf-8")


# -----------------------------
# Band extraction (symmetric, contiguous around ATM)
# -----------------------------
def band_from_df(df: pd.DataFrame, S: float):
    d = df.sort_values("Log-moneyness").reset_index(drop=True)
    i0 = (d["Log-moneyness"].abs()).idxmin()

    if not bool(d.loc[i0, "all_pass"]):
        return 0.0, (S, S)

    left = right = i0
    while left - 1 >= 0 and bool(d.loc[left - 1, "all_pass"]):
        left -= 1
    while right + 1 < len(d) and bool(d.loc[right + 1, "all_pass"]):
        right += 1

    x_left = float(d.loc[left, "Log-moneyness"])
    x_right = float(d.loc[right, "Log-moneyness"])
    x_star = min(abs(x_left), abs(x_right))

    K_min = S * math.exp(-x_star)
    K_max = S * math.exp(+x_star)
    return x_star, (K_min, K_max)


# -----------------------------
# Speed benchmark runner (creates speed_results.csv)
# -----------------------------
def run_speed_benchmark(S, r, T, sigma, out_csv, Ns=(100_000, 1_000_000), x_min=-0.1, x_max=0.1):
    coeff = atm_coeff(S, r, T, sigma)
    # handle 4 or 5 returns
    if len(coeff) == 4:
        C0, c1, c2, c3 = coeff
    else:
        C0, c1, c2, c3, _ = coeff

    rows = []
    for N in Ns:
        xs = np.linspace(x_min, x_max, num=int(N), dtype=np.float64)
        Ks = S * np.exp(-xs)

        t_bs = time_one(bs_vec, S, Ks, r, T, sigma)
        t_t3 = time_one(t3_vec, xs, S, C0, c1, c2, c3)
        speedup = t_bs / t_t3 if t_t3 > 0 else float("inf")

        rows.append({"N": int(N), "t_bs_sec": t_bs, "t_t3_sec": t_t3, "speedup": speedup})

    df = pd.DataFrame(rows)
    df.to_csv(out_csv, index=False)
    return df


# -----------------------------
# Main driver
# -----------------------------
def main():
    BASELINE_SWEEP_CSV = "outputs/baseline_sweep.csv"

    # ---- Update these if your filenames differ ----
    BASELINE_SWEEP_CSV = "baseline_sweep.csv"   # your sweep_strikes output saved to CSV
    GRID_BANDS_CSV     = "bands_grid.csv"       # optional grid output: columns [T, sigma, x_star]
    SPEED_CSV          = "speed_results.csv"    # will be created if missing

    os.makedirs("outputs", exist_ok=True)

    # ---- Baseline parameters (must match your baseline sweep file) ----
    S = 100.0
    r = 0.05
    T = 1.0
    sigma = 0.20

    # ---- Thresholds (set to match your project) ----
    PRICE_ABS = 0.01
    DELTA_ABS = 0.01
    GAMMA_REL = 0.05

    # ---- Load baseline sweep ----
    # ---- Baseline sweep: generate if missing ----
    if os.path.exists(BASELINE_SWEEP_CSV):
        baseline_df = pd.read_csv(BASELINE_SWEEP_CSV)
        if "gamma_rel_error" not in baseline_df.columns and "gamma_error" in baseline_df.columns:
            baseline_df = baseline_df.rename(columns={"gamma_error": "gamma_rel_error"})
    else:
        baseline_df = build_baseline_sweep_df(
            S, r, T, sigma,
            x_min=-0.05, x_max=0.05, step=0.001,
            save_path=BASELINE_SWEEP_CSV  # optional: saves it so future runs are instant
        )

    # Normalize gamma rel error column name if needed
    if "gamma_rel_error" not in baseline_df.columns and "gamma_error" in baseline_df.columns:
        baseline_df = baseline_df.rename(columns={"gamma_error": "gamma_rel_error"})

    # ---- Compute band ----
    x_star, (K_min, K_max) = band_from_df(baseline_df, S)

    # ---- Make baseline error plots ----
    plot_error_vs_x(
        baseline_df, x_star,
        ycol="call_error", ybound=PRICE_ABS,
        title="Call price absolute error vs log-moneyness",
        ylabel="|C_T3 - C_BS|",
        outpath="outputs/error_price.png"
    )

    plot_error_vs_x(
        baseline_df, x_star,
        ycol="delta_error", ybound=DELTA_ABS,
        title="Delta absolute error vs log-moneyness",
        ylabel="|Δ_T3 - Δ_BS|",
        outpath="outputs/error_delta.png"
    )

    plot_error_vs_x(
        baseline_df, x_star,
        ycol="gamma_rel_error", ybound=GAMMA_REL,
        title="Gamma relative error vs log-moneyness",
        ylabel="|Γ_T3 - Γ_BS| / |Γ_BS|",
        outpath="outputs/error_gamma_rel.png"
    )

    # ---- Heatmap (optional) ----
    if os.path.exists(GRID_BANDS_CSV):
        grid_df = pd.read_csv(GRID_BANDS_CSV)
        # requires columns: T, sigma, x_star
        if {"T", "sigma", "x_star"}.issubset(set(grid_df.columns)):
            heatmap_xstar(grid_df, outpath="outputs/heatmap_xstar.png")

    # ---- Speed results (generate if missing) ----
    if os.path.exists(SPEED_CSV):
        speed_df = pd.read_csv(SPEED_CSV)
    else:
        speed_df = run_speed_benchmark(S, r, T, sigma, out_csv=SPEED_CSV)

    # ---- Speed plots ----
    plot_speed(
        speed_df,
        outpath_times="outputs/speed_times.png",
        outpath_speedup="outputs/speed_speedup.png"
    )

    # ---- Summary txt ----
    params = {"S": S, "r": r, "T": T, "sigma": sigma}
    thresholds = {"price_abs": PRICE_ABS, "delta_abs": DELTA_ABS, "gamma_rel": GAMMA_REL}
    band_stats = {
        "x_star": x_star,
        "K_min": K_min,
        "K_max": K_max,
        "pct_minus": K_min / S - 1.0,
        "pct_plus":  K_max / S - 1.0,
    }
    speed_stats = speed_df.to_dict(orient="records")

    write_summary_txt("outputs/summary.txt", params, thresholds, band_stats, speed_stats)

    print("Done. Generated:")
    print("  outputs/error_price.png")
    print("  outputs/error_delta.png")
    print("  outputs/error_gamma_rel.png")
    if os.path.exists(GRID_BANDS_CSV):
        print("  outputs/heatmap_xstar.png (if grid file had correct columns)")
    print("  outputs/speed_times.png")
    print("  outputs/speed_speedup.png")
    print("  outputs/summary.txt")


if __name__ == "__main__":
    main()
