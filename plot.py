import matplotlib.pyplot as plt
import pandas as pd
import math
import numpy as np
from pathlib import Path


def plot_error_vs_x(df, x_star, ycol, ybound, title, ylabel, outpath):
    d = df.sort_values("Log-moneyness")
    x = d["Log-moneyness"].to_numpy()
    y = d[ycol].to_numpy()

    plt.figure()
    plt.plot(x, y)                 # don't set colors
    plt.axhline(ybound)            # threshold
    plt.axvline(-x_star)
    plt.axvline(+x_star)
    plt.title(title)
    plt.xlabel("log-moneyness x")
    plt.ylabel(ylabel)
    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    plt.close()

def band_summary_row(S, r, T, sigma, x_star, K_min, K_max, binder):
    return pd.DataFrame([{
        "S": S, "r": r, "T": T, "sigma": sigma,
        "x_star": x_star,
        "K_min": K_min, "K_max": K_max,
        "pct_minus": K_min/S - 1.0,
        "pct_plus":  K_max/S - 1.0,
        "binding_metric": binder
    }])

def heatmap_xstar(grid_df, outpath):
    pivot = grid_df.pivot(index="T", columns="sigma", values="x_star")
    plt.figure()
    plt.imshow(pivot.to_numpy(), aspect="auto", origin="lower")
    plt.title("Validated band radius x* across (T, sigma)")
    plt.xlabel("sigma grid index")
    plt.ylabel("T grid index")
    plt.xticks(range(len(pivot.columns)), [str(s) for s in pivot.columns])
    plt.yticks(range(len(pivot.index)), [str(t) for t in pivot.index])
    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    plt.close()

import matplotlib.pyplot as plt

def plot_speed(speed_df, outpath_times, outpath_speedup):
    d = speed_df.sort_values("N")

    # times
    plt.figure()
    plt.plot(d["N"], d["t_bs_sec"])
    plt.plot(d["N"], d["t_t3_sec"])
    plt.title("Wall-clock time vs batch size")
    plt.xlabel("N")
    plt.ylabel("seconds")
    plt.tight_layout()
    plt.savefig(outpath_times, dpi=200)
    plt.close()

    # speedup
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
    for k,v in params.items(): lines.append(f"  {k}: {v}\n")

    lines.append("\nThresholds:\n")
    for k,v in thresholds.items(): lines.append(f"  {k}: {v}\n")

    lines.append("\nValidated band (baseline):\n")
    for k,v in band_stats.items(): lines.append(f"  {k}: {v}\n")

    lines.append("\nSpeed benchmark:\n")
    for row in speed_stats:
        lines.append(f"  N={row['N']:,}: BS={row['t_bs_sec']:.4f}s, "
                     f"T3={row['t_t3_sec']:.4f}s, speedup={row['speedup']:.2f}x\n")

    p.write_text("".join(lines), encoding="utf-8")