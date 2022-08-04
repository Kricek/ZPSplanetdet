import numpy as np
from n_exp_maciek import get_massratio, get_N_exp
import emcee
import pandas as pd
import corner
import matplotlib.pyplot as plt
from tqdm import tqdm

"""
step by step analysis as described in Suzuki+16
pdf - probability density function
modyfikowana wersja mcmc_maciek.py z githuba
"""

def get_log_likelihood(params, q_break, S, sensitivity_data, planet_data):
    """
    pdf for the data given the parameters
    """
    log_s, log_q, surv_sens = sensitivity_data 
    s, q = planet_data
    A, m, n, p = params
    #nexp = get_N_exp(params, q_break, data)
    log_likelihood = -1 * get_N_exp(params, q_break, sensitivity_data) 
    for i in range(len(log_q)):    
        log_likelihood += np.log(A) + np.log(get_massratio(params, q_break, s[i], q[i]))
    return log_likelihood + np.log(S)

def get_log_prior(params):
    """
    uninformative prior pdf
    """
    A, m, n, p = params
    if  A > 0:
        return 0.0
    return -np.inf

def get_log_probability(params, q_break, S, sensitivity_data, planet_data):
    """
    posterior pdf for parameters given the data to be sampled
    """
    log_s, log_q, surv_sens = sensitivity_data
    s, q = planet_data
    lp = get_log_prior(params)
    if not np.isfinite(lp):
        return -np.inf
    return lp + get_log_likelihood(params, q_break, S, sensitivity_data, planet_data)
    
def calculate_S(data):
    log_s, log_q, surv_sens = data
    (max_j, max_i) = (surv_sens.shape) 
    S = 0
    for i in range(max_j):
        for j in range(max_j):    
            S += surv_sens[i][j]
    return S


# values imported from Suzuki+16

A = 0.61
m = 0.5
n = -0.92
p = 0.44
q_break = 1.7e-4
params = [A, m, n, p]

# reading data from survey_sensitivity2.dat

sensitivity_data = pd.read_csv('data/survey_sensitivity2.dat', header=None)
sensitivity_data = sensitivity_data.to_numpy()
log_s = sensitivity_data[0,2:]
log_q = sensitivity_data[1:,1]
surv_sens = sensitivity_data[1:,2:]
sensitivity_data = [log_s, log_q, surv_sens]

# reading data from planets_data.dat

planet_data = pd.read_csv('data/planets_data.csv', header=None)
planet_data = planet_data.to_numpy()
s = planet_data[1:,2].astype(float)
q = planet_data[1:,1].astype(float)
q *= 10**-3
planet_data = [s, q]

# calculating S

S = calculate_S(sensitivity_data)

# setting starting positions for walkers

pos = np.random.rand(15, 4)
nwalkers, ndim = pos.shape

# running emcee

sampler = emcee.EnsembleSampler(nwalkers, ndim, get_log_probability, args = (q_break, S, sensitivity_data, planet_data))
sampler.run_mcmc(pos, 1000, progress=True)

# visualization of results

labels = ["A", "m", "n", "p"]
flat_samples = sampler.get_chain(flat = True)
figure = corner.corner(flat_samples, labels = labels)
plt.show()
