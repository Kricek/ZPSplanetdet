import numpy as np
import pandas as pd

def get_massratio(s, q, q_break, params): 
    """
    generating mass-ratio function in shape of f(s, q, params = [q_break,n,p,m]) = [(q/q_break)^n+(q/q_break)^p]*s^m 
    params = [q_break, n, p, m]
    """
    n, p, m = params
    if q > q_break:        
        return (q/q_break)**n*s**m
    else:
        return (q/q_break)**p*s**m

def get_N_exp(params, q_break, data, Aqr): 
    """
    Trivial integral by summation through all data given of function f*S where f is mass-ration function, S is survey sensitivity given by data
    params = [A, q_break, n, p, m]
    data = [log_s, log_q, survey_sensitivity]
    Aqr = A constant in mass ratio distribution
    """
    integral = 0
    log_s, log_q, surv_sens = data
    for i in range(len(log_q)-1):
        for j in range(len(log_s)-1):
            dsdq = (log_s[j+1]-log_s[j])*(log_q[i+1]-log_q[i])
            integral += dsdq*get_massratio(10**log_s[j],10**log_q[i], q_break, params)*surv_sens[i,j]
    return Aqr*integral

sensitivity_data = pd.read_csv('data/survey_sensitivity2.dat', header=None)
data = sensitivity_data.to_numpy()
log_s = data[0,2:]
log_q = data[1:,1]
surv_sens = data[1:,2:]
A = 0.62
q_break = 1.7e-4
n = -0.92
p = 0.47
m = 0.5
params = [n, p, m]
N_exp = get_N_exp(params, q_break, [log_s, log_q, surv_sens], A)
print(N_exp)
