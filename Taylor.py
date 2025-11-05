import math
import numpy as np
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

price_error_bound = 0.01
price_bps_error_bound = 2.0
delta_error_bound = 0.005
gamma_error_bound = 0.005  #

def errors_at_x(s,r,t,sigma,x,c0,c1,c2,c3):
    k = strike(s,x)
    call_exact = bsm_call_value(s, k, r, t, sigma)
    delta_exact = bsm_delta(s, k, r, t, sigma)
    gamma_exact = bsm_gamma(s, k , r ,t, sigma)

    call_taylor = taylor_price_x(x, c0, c1, c2, c3)
    delta_taylor = delta_price_x(x, s, c1, c2, c3)
    gamma_taylor = gamma_price_x(x, s, c1, c2, c3)

    call_error = abs(call_exact - call_taylor)
    delta_error = abs(delta_exact - delta_taylor)
    gamma_error = abs(gamma_exact - gamma_taylor)

    price_pass = call_error <= price_error_bound
    delta_pass = delta_error <= delta_error_bound
    gamma_pass = gamma_error <= gamma_error_bound

    return { "x" : x, "k" : k, "call_exact" : call_exact,
              "delta_exact" : delta_exact, "gamma_exact" : gamma_exact,
              "call_taylor" : call_taylor, "delta_taylor" : delta_taylor,
              "gamma_taylor" : gamma_taylor, "call_error" : call_error,
              "delta_error" : delta_error, "gamma_error" : gamma_error,
             "pass_all" : (price_pass and delta_pass and gamma_pass)}

#class to sweep multiple x values (different strikes) at makes a list
def sweep_strikes(s,r,t,sigma, x_min, x_max, step):
    c0, c1, c2, c3, c4 = atm_coeff(s,r,t,sigma)
    strikes = list(np.arange(x_min, x_max + step, step))
    taylor_price_list = []
    for x in strikes:
        taylor_price_list.append(taylor_price_x(x, c0, c1, c2, c3))

    return taylor_price_list


coeff = atm_coeff(s, r, t, sigma)
c0 = coeff[0]
c1 = coeff[1]
c2 = coeff[2]
c3 = coeff[3]
c4 = coeff[4]

k = strike(s, x)

call_value = taylor_price_x(x, c0, c1, c2, c3)
call_value_bsm = bsm_call_value(s, k, r, t, sigma)
#print("x = " + str(x))
#print("Taylor value: " + str(call_value))
#print("Bsm value: " + str(call_value_bsm))
#print(c4)

print(sweep_strikes(s,r,t,sigma, -0.1,0.1,0.05))
print(bsm_call_value(s,105.127,r,t,sigma))