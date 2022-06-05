import numpy as np
from n_exp_maciek import get_N_exp, get_massratio, get_sensitivity
import emcee
import pandas as pd
import corner
import matplotlib.pyplot as plt
from tqdm import tqdm

"""
step by step analysis as described in Suzuki+16
pdf - probability density function
"""

def get_log_likelihood(params, q_break, data):
    """
    pdf for the data given the parameters
    """
    s, q, surv_sens = data # s and q are decimal logarithms
    log_likelihood = get_N_exp(params, q_break, data)
    for i in range(len(q)):
        log_likelihood += np.log(get_massratio(s[i], q[i], q_break, params)) + np.log(get_sensitivity(s[i], q[i], data))
    return log_likelihood

def get_log_prior(params):
    """
    uninformative prior pdf
    """
    A, n, p, m = params
    if  A > 0:
        return 0.0
    return -np.inf

def get_log_probability(params, q_break, data):
    """
    posterior pdf for parameters given the data to be sampled
    """
    s, q, surv_sens = data
    lp = get_log_prior(params)
    if not np.isfinite(lp):
        return -np.inf
    return lp + get_log_likelihood(params, q_break, data)


#values imported from Suzuki+16

A = 0.62
q_break = 1.7e-4
n = -0.92
p = 0.47
m = 0.5
params = [A, n, p, m]

#reading data from survey_sensitivity2.dat

sensitivity_data = pd.read_csv('data/survey_sensitivity2.dat', header=None)
data = sensitivity_data.to_numpy()
s = data[0,2:]
q = data[1:,1]
surv_sens = data[1:,2:]
surv_data = [s, q, surv_sens]

#setting starting positions and number of walkers

pos = np.random.rand(15, 4)
nwalkers, ndim = pos.shape

#running emcee

sampler = emcee.EnsembleSampler(nwalkers, ndim, get_log_probability, args = (q_break, surv_data))
sampler.run_mcmc(pos, 1000, progress=True)

#visualization of results

flat_samples = sampler.get_chain(flat = True)
figure = corner.corner(flat_samples)
plt.show()
