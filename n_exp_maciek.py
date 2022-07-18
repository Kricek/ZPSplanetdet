import numpy as np
import pandas as pd

def get_massratio(params, q_break, s, q): 
    """
    generating mass-ratio function in shape of f(s, q, params = [A,q_break,m,n,p]) = A*[(q/q_break)^n+(q/q_break)^p]*s^m 
    params = [A, m, n, p]
    """
    A, m, n, p = params
    if q > q_break:   
        result = ((q/q_break)**n)*(s**m)
        #print(result)      
        return result
    else:
        result = ((q/q_break)**p)*(s**m)
        #print(result)      
        return result

def get_N_exp(params, q_break, data): 
    """
    trivial integral by summation through all data given of function f*S where f is mass-ration function, S is survey sensitivity given by data
    params = [A, q_break, m, n, p]
    data = [log_s, log_q, survey_sensitivity]
    A = constant in mass ratio distribution
    """
    integral = 0
    log_s, log_q, surv_sens = data
    A, m, n, p = params
    for i in range(len(log_q)-1):
        for j in range(len(log_s)-1):
            dsdq = (log_s[j+1]-log_s[j])*(log_q[i+1]-log_q[i])
            integral += dsdq*get_massratio(params, q_break, 10**log_s[j],10**log_q[i])*surv_sens[i,j]
    return A * integral

sensitivity_data = pd.read_csv('data/survey_sensitivity2.dat', header=None)
data = sensitivity_data.to_numpy()
log_s = data[0,2:]
log_q = data[1:,1]
surv_sens = data[1:,2:]

A = 0.62
m = 0.5
n = -0.92
p = 0.47
q_break = 1.7e-4
params = [A, m, n, p]

N_exp = get_N_exp(params, q_break, [log_s, log_q, surv_sens])
print(N_exp)
