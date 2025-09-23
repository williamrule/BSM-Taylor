import math
import numpy as np
from scipy.stats import norm
from BSM import bsm_call_value

#Constant paramaters for testing
s = 100
r = 0.05
t = 1
sigma = 0.20

def coeff_bounds(s, r, t, sigma, x):
    #note, we are atm so s = k so no need to call strike
    a = sigma * math.sqrt(t)
    m = x + (r + 0.5 * sigma ** 2) * t
    d1 =  m / a
    d2 = d1 - a

    delta = norm.cdf(d1)

    c0 = (s * delta) - (s * np.exp(-r * t) * norm.cdf(d2))
    c1 = s * delta
    c2 = c1 + (s * norm.pdf(d1))/a
    c3 = c1 + (2 * s * norm.pdf(d1))/a - (s * d1 * norm.pdf(d1))/ a ** 2

    return c0, c1, c2, c3

coeff = coeff_bounds(s, r, t, sigma, 0)
c0 = coeff[0]
c1 = coeff[1]
c2 = coeff[2]
c3 = coeff[3]

print(coeff[3])