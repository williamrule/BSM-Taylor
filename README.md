# BSM-Taylor
**Taylor-series approximation of the Black–Scholes–Merton (BSM) European option price**, with error + speed experiments to map *where the approximation holds* and *where it breaks*.

> **Core idea:** Expand the BSM price as a Taylor series around an expansion point (often near ATM), then quantify approximation error as you move away in **moneyness**, **volatility**, and **time to maturity**—and compare runtime vs. the closed-form pricer.

---

## Why this project
Closed-form BSM is already fast, but in real quant workflows you often care about:
- **Approximation behavior** (stability, error surfaces, breakdown regimes)
- **Fast repeated evaluation** (calibration loops, scenario sweeps, Greeks/proxies)
- **Understanding model sensitivity** via controlled perturbations

This repo is a compact, experiment-driven exploration of those themes.

---

## What’s in this repo (current layout)
- `BSM.py` — Black–Scholes(-Merton) pricing functions (baseline “ground truth”)
- `Taylor.py` — Taylor approximation logic built around an expansion point
- `error.py` — error metrics (absolute/relative), comparisons vs. exact BSM
- `sweep.py` — parameter sweeps over spot/moneyness/σ/T (generates results)
- `speed.py` — timing comparisons (exact BSM vs. Taylor approximation)

> Tip: If you’re viewing this as a portfolio repo, the README + plots matter more than extra code.

---

## Math (high level)
For a European call, the BSM price is:

\[
C(S,K,r,\sigma,T) = S \, N(d_1) - K e^{-rT} N(d_2)
\]
\[
d_1 = \frac{\ln(S/K) + (r + \tfrac{1}{2}\sigma^2)T}{\sigma \sqrt{T}}, \quad
d_2 = d_1 - \sigma \sqrt{T}
\]

A Taylor approximation expands the price around an expansion point \(S_0\) (or another variable such as log-moneyness), e.g.

\[
C(S) \approx \sum_{n=0}^{N} \frac{C^{(n)}(S_0)}{n!}(S - S_0)^n
\]

This repo focuses on **empirical error behavior** as you move away from the expansion point and vary parameters.

---

## Setup
### Dependencies
You’ll typically need:
- `numpy`
- `scipy` (if you use normal CDF/PDF from SciPy)
- `matplotlib` (for plots)

Install quickly:
```bash
pip install numpy scipy matplotlib
```

(If you later add a `requirements.txt`, replace the above with `pip install -r requirements.txt`.)

---

## Quickstart
Run the scripts directly:
```bash
python sweep.py
python error.py
python speed.py
```

If scripts require editing parameters, look near the top of each file for constants like:
- strike `K`
- rate `r`
- volatility `sigma`
- maturity `T`
- expansion point `S0`
- polynomial order `N`
- sweep grid ranges / step sizes

---

## Experiments you can run
### 1) Error vs. spot (moneyness sweep)
Goal: Visualize how error grows as \(S\) moves away from the expansion point.

What to look for:
- Error is typically smallest near \(S_0\) and grows with distance
- Higher Taylor order often improves the local region but can still blow up further out

### 2) Error surface / heatmap
Goal: Sweep across two dimensions (e.g., moneyness × volatility).

What to look for:
- Regions where approximation is reliable (low error “basin”)
- Parameter regimes where it breaks sharply (high vol, long maturities, deep ITM/OTM)

### 3) Speed vs. accuracy tradeoff
Goal: Compare runtime of exact BSM vs. Taylor approximation across many evaluations.

What to look for:
- Whether the approximation provides speed gains at tolerable error
- How runtime scales with Taylor order \(N\)

---

## Recommended outputs (for a portfolio-ready repo)
Add a `figures/` folder and include 2–4 of your best plots here. For example:
- `figures/error_vs_spot.png`
- `figures/error_heatmap.png`
- `figures/speed_comparison.png`

Then embed them:

### Error vs. Spot
![Error vs Spot](figures/error_vs_spot.png)

### Error Heatmap (Moneyness × Volatility)
![Error Heatmap](figures/error_heatmap.png)

### Speed Comparison
![Speed Comparison](figures/speed_comparison.png)

> If your scripts don’t currently save figures automatically: add a `plt.savefig("figures/<name>.png", dpi=200, bbox_inches="tight")`.

---

## Interpretation notes (what this repo is trying to show)
- Taylor approximations are **local**: accuracy is best near the expansion point.
- Approximation quality depends strongly on **moneyness**, **σ**, and **T**.
- Increasing polynomial order can help locally but may be unstable far away.
- The interesting part is not “Taylor works” — it’s **mapping the boundary** of where it stops working and why.

---

## Limitations / scope
- Educational / research prototype quality (not production-ready pricing code)
- Assumes European options under standard BSM assumptions
- Does not attempt to model jumps, stochastic volatility, dividends (unless explicitly implemented)

---

## Roadmap (nice upgrades if you want)
- Package structure (`src/bsm_taylor/…`) and importable modules
- CLI arguments (e.g., `python sweep.py --K 100 --sigma 0.2 --T 1.0 --order 6`)
- Unit tests (sanity: monotonicity in S, call-put parity, boundary behavior)
- Adaptive order selection (choose \(N\) based on tolerance + distance from \(S_0\))
- Expand around **log-moneyness** for improved numerical behavior

---

## References
- Black, F. & Scholes, M. (1973). *The Pricing of Options and Corporate Liabilities.*
- Merton, R. C. (1973). *Theory of Rational Option Pricing.*

---

## License
Add a license if you want others to reuse this (MIT is common for portfolios).
