# BSM-Taylor
Taylor-series approximation of the **Black–Scholes–Merton (BSM)** European option price, plus experiments that measure **accuracy (price/Δ/Γ)** and **speed** vs. the exact closed-form BSM pricer.

**What this repo demonstrates:** approximation behavior is *local*—it can be very accurate near the expansion point (often near-ATM), and error grows as you move away in moneyness / volatility / maturity.

---

## What’s included

### Core modules
- `BSM.py` — baseline BSM pricing / Greeks used as the “ground truth”
- `Taylor.py` — Taylor approximation implementation around an expansion point

### Experiments / utilities
- `sweep.py` — parameter sweep(s) to compare Taylor vs exact BSM
- `error.py` — error computation (price / delta / gamma)
- `speed.py` — timing comparisons (exact vs Taylor)
- `plot.py` — plotting utilities for the experiments
- `make_outputs.py` — convenience script to regenerate the artifacts in `outputs/` (CSV + plots)

---

## Outputs (already in this repo)

The `outputs/` directory contains:
- `baseline_sweep.csv` — the sweep results used to produce the plots
- `error_price.png` — pricing error plot
- `error_delta.png` — delta error plot
- `error_gamma_rel.png` — relative gamma error plot
- `speed_times.png` — timing comparison plot
- `speed_speedup.png` — speedup plot

> Note: `outputs/summary.txt` is a run log / summary and is meant to be ignored (not tracked).

### Preview

#### Price error
![Price error](outputs/error_price.png)

#### Delta error
![Delta error](outputs/error_delta.png)

#### Gamma error (relative)
![Gamma error (relative)](outputs/error_gamma_rel.png)

#### Timing
![Timing](outputs/speed_times.png)

#### Speedup
![Speedup](outputs/speed_speedup.png)

---

## Math (high level)

BSM European call price:

$$
C(S,K,r,\sigma,T) = S\,N(d_1) - K e^{-rT} N(d_2)
$$

$$
d_1 = \frac{\ln(S/K) + (r + \tfrac{1}{2}\sigma^2)T}{\sigma\sqrt{T}},\quad
d_2 = d_1 - \sigma\sqrt{T}
$$

Taylor approximation idea (illustrative form, expanding in spot around \(S_0\)):

$$
C(S) \approx \sum_{n=0}^{N}\frac{C^{(n)}(S_0)}{n!}(S-S_0)^n
$$

This repo focuses on **empirical error behavior** (how approximation quality changes across parameter regimes) and **runtime tradeoffs** (speed vs. accuracy).

---

## Setup

### Dependencies
This project uses standard scientific Python tooling. Install:
```bash
pip install numpy scipy matplotlib
