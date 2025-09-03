import math
import numpy as np
from scipy.stats import norm

#s = Stock Price, k = Strike Price, r = Risk free interest rate, sigma = volatility

#Exact BSM formula
def bsm_call_value(s, k, r, t,sigma):
    d1 = (math.log(s / k) + (r + 0.5 * sigma**2) * t) / sigma * math.sqrt(t)
    d2 = d1 - sigma * math.sqrt(t)

    c = (s * norm.cdf(d1)) - (k* np.exp(-r * t) * norm.cdf(d2))

    return c

#Delta: First derivitive of c with respect to stock price
def bsm_delta(s, k, r, t, sigma):
    d1 = (math.log(s / k) + (r + 0.5 * sigma ** 2) * t) / sigma * math.sqrt(t)
    return s * norm.cdf(d1)

def bsm_gamma(s, k, r, t, sigma):
    d1 = (math.log(s / k) + (r + 0.5 * sigma ** 2) * t) / sigma * math.sqrt(t)
    first = s * norm.cdf(d1)
    second = s * norm.pdf(first) / (sigma * math.sqrt(t))
    return first + second

print(bsm_call_value(100, 100, 0.05, 1.0, 0.20))
print(bsm_delta(100, 100, 0.05, 1.0, 0.20))
print(bsm_gamma(100, 100, 0.05, 1.0, 0.20))

#Need to fix gamma, wrong formula