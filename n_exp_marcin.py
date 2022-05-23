import pandas as pd
import numpy as np

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
    
def get_N_exp(s, q, sensitivity, q_break, Aqr):
    dV = np.zeros((len(log_q),len(log_s)))
    for i in range(len(log_q)):
        if i == 0:
            dq = (q[i+1]-q[i])
        elif i == len(log_q) - 1:
            dq = (q[i]-q[i-1])
        else:
            dq = (q[i]-q[i-1])/2 + (q[i+1]-q[i])/2
                     
        for j in range(len(log_s)):
            if j == 0:
                ds = (s[j+1]-s[j])
            elif j == len(log_s) - 1:
                ds = (s[j]-s[j-1])
            else:
                ds = (s[j]-s[j-1])/2 + (s[i+1]-s[i])/2
            dV[i,j] = ds * dq * surv_sens[i,j] * get_massratio(s[j], q[i], q_break, params)
    return Aqr*np.sum(dV)
    
    
sensitivity_data = pd.read_csv('data/survey_sensitivity2.dat', header=None)
data = sensitivity_data.to_numpy()
log_s = data[0,2:]
log_q = data[1:,1]
surv_sens = data[1:,2:]

s = 10**log_s
q = 10**log_q
A = 0.62
q_break = 1.7e-4
n = -0.92
p = 0.47
m = 0.5
params = [n, p, m]

N_exp = get_N_exp(s, q, surv_sens, q_break, A)
print(N_exp)
    
