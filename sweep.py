import numpy as np
import pandas as pd
from Taylor import atm_coeff
from error import errors_at_x

#Constant parameters for testing
s = 100
r = 0.05
t = 1
sigma = 0.20


#class to sweep multiple x values (different strikes) and makes a list with the taylor price and error at that price
def sweep_strikes(s,r,t,sigma, x_min, x_max, step):
    c0, c1, c2, c3, c4 = atm_coeff(s,r,t,sigma)
    strikes = list(np.arange(x_min, x_max + step, step))
    taylor_strike_list = []
    for i in strikes:
        taylor_strike_list.append(errors_at_x(s,r,t,sigma, i, c0, c1, c2, c3))
    return taylor_strike_list

#Returns the band that satisfies the desired error and a dataframe of the passed strikes
#Note to self: Need to make band symmetric and continguity

def band_extraction(taylor_strike_list):
    df = pd.DataFrame(taylor_strike_list)
    all_pass_list = df[df["all_pass"]]
    col_min = all_pass_list["Log-moneyness"].min()
    col_max = all_pass_list["Log-moneyness"].max()
    band = [col_min, col_max]
    return all_pass_list, band

taylor_strike_list = sweep_strikes(s,r,t,sigma, -0.05,0.050,.001)
all_pass_list, band = band_extraction(taylor_strike_list)
print(all_pass_list.info())
print(all_pass_list)
print(band)


'''
#Sweep test
strike_list = sweep_strikes(s,r,t,sigma, -0.05,0.050,.001)
df = pd.DataFrame(strike_list)
print(df.info())
print(df[df['all_pass'] == True])
df.to_csv("new_test")
'''