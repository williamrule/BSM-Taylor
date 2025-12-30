from scipy.stats import norm
from BSM import bsm_call_value, bsm_delta, bsm_gamma
from Taylor import taylor_price_x, strike, delta_price_x, gamma_price_x

#Constant parameters for testing
s = 100
r = 0.05
t = 1
sigma = 0.20

#pre-determined error bounds
price_error_bound = 0.1
price_bps_error_bound = 2.0
delta_error_bound = 0.01
gamma_rel_bound = 0.05

#computes errors at x from Taylor and BSM and makes a dict
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
    gamma_rel_error = abs(gamma_exact - gamma_taylor) / max(abs(gamma_exact), 1e-16)


    price_pass = call_error <= price_error_bound
    delta_pass = delta_error <= delta_error_bound
    gamma_pass = gamma_rel_error <= gamma_rel_bound

    return {"Log-moneyness" : x,
            "strike" : k,
            "call_exact" : call_exact,
            "call_taylor": call_taylor,
            "call_error": call_error,
            "delta_exact": delta_exact,
            "delta_taylor": delta_taylor,
            "delta_error": delta_error,
            "gamma_exact": gamma_exact,
            "gamma_taylor": gamma_taylor,
            "gamma_error": gamma_rel_error,
            "call_pass": price_pass,
            "all_pass": (price_pass and delta_pass and gamma_pass)
            }
