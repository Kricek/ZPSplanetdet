import numpy as np
import scipy.integrate as integrate

def fun_massratio(s, q, params): 
    """
    generating mass-ratio function in shape of f(s,q,A,n,m) = A*q^n*s^m 
    params = [A, n, m]
    """
    return params[0]*q**params[1]*s**params[2]

def fun_sensitivity(s,q): #should be moved to read_survey_sensitivity file as function interpolating S(s,q)
    return .2

def fun_massratio_sensitivity(s, q, params): #function to integrate: f*S
    return fun_massratio(s,q, params)*fun_sensitivity(s,q)

def get_N_exp(params): #integration with finite limits q in (0,1), s in (0,5)
    return integrate.dblquad(fun_massratio_sensitivity, 0, 1, 0, 5, args=(params,))

A = 1
n = 2
m = 2
params = [A, n, m]
N_exp, N_exp_error = get_N_exp(params)
print(N_exp)