import numpy as np
from n_exp import get_N_exp, get_massratio, get_sensitivity
import emcee
import pandas as pd

"""
step by step analysis as described in Suzuki+16
pdf - probability density function
"""

def get_log_likelihood(params, survey_data, planet_data):
    """
    pdf for the data given the parameters
    """
    s_i, q_i = planet_data
    log_likelihood = get_N_exp(params, survey_data)
    for i in range(len(s_i)):
        log_likelihood *= np.log(get_massratio(s_i[i], q_i[i], params)) + np.log(get_sensitivity(s_i[i], q_i[i], survey_data))
    return log_likelihood

def get_log_prior(params):
    """
    uninformative prior pdf
    """
    A, n, p, m = params
    if  A > 0:
        return 0.0
    return -np.inf

def get_log_probability(params, survey_data, planet_data):
    """
    posterior pdf for parameters given the data to be sampled
    """
    lp = get_log_prior(params)
    if not np.isfinite(lp):
        return -np.inf
    return lp + get_log_likelihood(params, survey_data, planet_data)

params = (A, n, p, m)
pos = np.random(size=(15, 4))
nwalkers, ndim = pos.shape

sampler = emcee.EnsembleSampler(
nwalkers, ndim, get_log_probability, args = ()
)
sampler.run_mcmc(pos, 1000, progress=True)















