import math
import numpy as np
import pandas as pd
from scipy.stats import norm
from BSM import bsm_call_value, bsm_delta, bsm_gamma

#Constant paramaters for testing
s = 100
r = 0.05
t = 1
sigma = 0.20
x = 0.15

#Define everything in one call
def atm_coeff(s, r, t, sigma):
    #note, we are atm so s = k so no need to call strike
    a = sigma * math.sqrt(t)
    m = (r + 0.5 * sigma ** 2) * t
    d1 =  m / a
    d2 = d1 - a

    delta = norm.cdf(d1)

    c0 = (s * delta) - (s * np.exp(-r * t) * norm.cdf(d2))
    c1 = s * delta
    c2 = c1 + (s * norm.pdf(d1))/a
    c3 = c1 + (2 * s * norm.pdf(d1))/a - (s * d1 * norm.pdf(d1))/ a ** 2
    c4 = c1 + (3 * s * norm.pdf(d1))/a - (3 * d1 * norm.pdf(d1))/ a ** 2 - (norm.pdf(d1)) / a ** 3 + (d1 ** 2) + norm.pdf(d1) / a ** 3

    return c0, c1, c2, c3, c4

#Compute taylor price from x
#Remember x = ln(s/k)
def taylor_price_x(x, c0, c1, c2, c3):
    return c0 + (x * (c1 + (x * (1/2 * c2 + (x * 1/6 * c3)))))

def delta_price_x(x, s, c1, c2, c3):
    return (c1 + (c2 * x) + (0.5 * c3 * (x**2)))/ s
#Gamma approximation
def gamma_price_x(x, s, c1, c2, c3):
    return ((c2 - c1) + ((c3 - c2) * x)  - (0.5 * c3 * (x**2)))/ (s**2)

def strike(s, x):
    return s * math.exp(-x)




coeff = atm_coeff(s, r, t, sigma)
c0 = coeff[0]
c1 = coeff[1]
c2 = coeff[2]
c3 = coeff[3]
c4 = coeff[4]



