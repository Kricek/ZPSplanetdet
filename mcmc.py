import numpy as np
from n_exp import get_N_exp, fun_massratio, fun_sensitivity
import emcee
import pandas as pd

def get_log_likelihood(params, survey_data, planet_data):
    s_i, q_i = planet_data
    likelihood =  np.exp(get_N_exp(params, survey_data))
    for i in range(len(s_i)):
        likelihood *= fun_massratio(s_i[i], q_i[i], params)*fun_sensitivity(s_i[i], q_i[i], survey_data)
    return np.log10(likelihood)

def get_log_prior(params):
    A, q_break, n, p, m = params
    if  A > 0 and q_break > 0 and q_break < 1:
        return 0.0
    return -np.inf

def log_probability(params, survey_data, planet_data):
    lp = get_log_prior(params)
    if not np.isfinite(lp):
        return -np.inf
    return lp + get_log_likelihood(params, survey_data, planet_data)

