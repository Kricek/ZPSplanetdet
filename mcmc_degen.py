import numpy as np
from n_exp_maciek import get_massratio, get_N_exp
from surv_ip_kuba import InterpolateSens
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

def get_log_likelihood(params, q_break, sensitivity_data, planet_data):
    """
    pdf for the data given the parameters
    """
    log_s, log_q, surv_sens = sensitivity_data 
    s, q, w = planet_data
    A, m, n, p = params
    log_L1, log_L2 = (0, 0)
    #nexp = get_N_exp(params, q_break, data)
    log_likelihood = -1 * get_N_exp(params, q_break, sensitivity_data) 
    for i in range(len(q)):
        if(w[i] == 1):    
            log_L1 += np.log(A) + np.log(get_massratio(params, q_break, s[i], q[i])) + np.log(InterpolateSens(log_s, log_q, surv_sens, s[i], q[i], 4, 'cubic'))
        else:
            L2a = A * get_massratio(params, q_break, s[i], q[i]) * InterpolateSens(log_s, log_q, surv_sens, s[i], q[i], 4, 'cubic') * w[i]
            L2b = A * get_massratio(params, q_break, s[i+1], q[i+1]) * InterpolateSens(log_s, log_q, surv_sens, s[i+1], q[i+1], 4, 'cubic') * w[i+1]
            log_L2 += np.log(L2a + L2b)
            i += 2
        if(i > 25):
            break 
    log_likelihood += log_L1 + log_L2 
    return log_likelihood # + np.log(S)

def get_log_prior(params):
    """
    uninformative prior pdf
    """
    A, m, n, p = params
    if  A > 0:
        return 0.0
    return -np.inf

def get_log_probability(params, q_break, sensitivity_data, planet_data):
    """
    posterior pdf for parameters given the data to be sampled
    """
    log_s, log_q, surv_sens = sensitivity_data
    s, q, w = planet_data
    lp = get_log_prior(params)
    if not np.isfinite(lp):
        return -np.inf
    return lp + get_log_likelihood(params, q_break, sensitivity_data, planet_data)

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
w = planet_data[1:,3].astype(float)
q *= 10**-3
planet_data = [s, q, w]

# setting starting positions for walkers

pos = np.random.rand(50, 4)
nwalkers, ndim = pos.shape

# running emcee

sampler = emcee.EnsembleSampler(nwalkers, ndim, get_log_probability, args = (q_break, sensitivity_data, planet_data))
sampler.run_mcmc(pos, 10000, progress=True)

# visualization of results

labels = ["A", "m", "n", "p"]
flat_samples = sampler.get_chain(flat = True)
figure = corner.corner(flat_samples, labels = labels)
plt.show()
Footer

