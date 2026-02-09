import numpy as np
from scipy.stats import norm
import time, gc, statistics
from Taylor import atm_coeff

# Baseline Parameter set
S = 100.0
r = 0.05
T = 1.0
sigma = 0.20

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

C0, c1, c2, c3, c4 = atm_coeff(S, r, T, sigma)

for N in [100_000, 1_000_000]:
    xs = np.linspace(-0.1, 0.1, num=N)
    Ks = S * np.exp(-xs)

    t_bs = time_one(bs_vec, S, Ks, r, T, sigma)
    t_t3 = time_one(t3_vec, xs, S, C0, c1, c2, c3)
    speedup = t_bs / t_t3

    print(f"N={N:,} | BS={t_bs:.4f}s | T3={t_t3:.4f}s | speedup={speedup:.2f}x")
