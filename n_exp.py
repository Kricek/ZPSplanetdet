import numpy as np
import pandas as pd

def get_lower_nearest_index(array, x): 
    """
    giving index of array element closest to x
    """
    idx = (np.abs(array - x)).argmin()
    if array[idx] > x:
        return idx - 1
    else:
        return idx

def fun_massratio(s, q, q_break, params): 
    """
    generating mass-ratio function in shape of f(s, q, params = [A,q_break,n,p,m]) = A*[(q/q_break)^n+(q/q_break)^p]*s^m 
    params = [q_break, n, p, m]
    """
    n, p, m = params
    if q > q_break:        
        return (q/q_break)**n*s**m
    else:
        return (q/q_break)**p*s**m

def fun_sensitivity(s,q, data):
    """
    interpolate survey sensitivity from data given as arrays (log_s, log_q, surv_sens) for given point (s,q)
    """
    log_s, log_q, surv_sens = data
    s_low_idx = get_lower_nearest_index(log_s, np.log10(s))
    q_low_idx = get_lower_nearest_index(log_q, np.log10(q))

    if s_low_idx == 0 or q_low_idx == 0 or s_low_idx == (len(log_s)-1) or q_low_idx == (len(log_q)-1):
        interpolated_surv_sens = surv_sens[q_low_idx,s_low_idx]
        return interpolated_surv_sens
    
    point_A = [log_q[q_low_idx], log_s[s_low_idx], surv_sens[q_low_idx, s_low_idx]]
    point_B = [log_q[q_low_idx+1], log_s[s_low_idx], surv_sens[q_low_idx+1, s_low_idx]]
    point_C = [log_q[q_low_idx], log_s[s_low_idx+1], surv_sens[q_low_idx, s_low_idx+1]]
#    print(point_A)
#    print(point_B)
#    print(point_C)
    vector_AB = np.array(point_B) - np.array(point_A)
    vector_AC = np.array(point_C) - np.array(point_A)
    vector_normal = np.cross(vector_AB, vector_AC)
    plane_constant = np.dot(vector_normal, point_A)
    interpolated_surv_sens = (-vector_normal[0]*np.log10(q)-vector_normal[1]*np.log10(s)+plane_constant)/vector_normal[2]
#    print([np.log10(q), np.log10(s), interpolated_surv_sens])
    return interpolated_surv_sens

def get_N_exp(params, q_break, data, Aqr): 
    """
    Trivial integral by summation through all data given of function f*S where f is mass-ration function, S is survey sensitivity given by data
    params = [A, q_break, n, p, m]
    data = [log_s, log_q, survey_sensitivity]
    Aqr = A constant in 
    """
    integral = 0
    log_s, log_q, surv_sens = data
    for i in range(len(log_q)-1):
        for j in range(len(log_s)-1):
            dsdq = (log_s[i+1]-log_s[i])*(log_q[i+1]-log_q[i])
            integral += dsdq*fun_massratio(10**log_s[i],10**log_q[i], q_break, params)*surv_sens[i,j]
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