import numpy as np
import emcee
import pandas as pd
import corner
import matplotlib.pyplot as plt
from tqdm import tqdm
from datetime import datetime

from n_exp_maciek import get_massratio, get_N_exp
from surv_ip_kuba import InterpolateSens

"""
step by step analysis as described in Suzuki+16
pdf - probability density function
"""

def get_log_likelihood(params, q_break, sensitivity_data, planet_data, S):
    """
    pdf for the data given the parameters
    """
    log_s, log_q, surv_sens = sensitivity_data 
    s, q, w = planet_data
    A, m, n, p = params
    log_L1, log_L2 = (0, 0)
    #nexp = get_N_exp(params, q_break, data)
    log_likelihood = -1 * get_N_exp(params, q_break, sensitivity_data)
    i = 0
    while i < len(q):        
        if(w[i] == 1):    
            log_L1 += np.log(A) + np.log(get_massratio(params, q_break, s[i], q[i])) + np.log(S[i])
            i += 1
        else:
            L2a = A * get_massratio(params, q_break, s[i], q[i]) * S[i] * w[i]
            L2b = A * get_massratio(params, q_break, s[i+1], q[i+1]) * S[i+1] * w[i+1]
            log_L2 += np.log(L2a + L2b)
            i += 2
    log_likelihood += log_L1 + log_L2 
    return log_likelihood 

def get_log_prior(params):
    """
    uninformative prior pdf
    """
    A, m, n, p = params
    if  A > 0:
        return 0.0
    return -np.inf

def get_log_probability(params, q_break, sensitivity_data, planet_data, S):
    """
    posterior pdf for parameters given the data to be sampled
    """
    log_s, log_q, surv_sens = sensitivity_data
    s, q, w = planet_data
    lp = get_log_prior(params)
    if not np.isfinite(lp):
        return -np.inf
    return lp + get_log_likelihood(params, q_break, sensitivity_data, planet_data, S)

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

# calculating S for all (s_i, q_i)

S = np.zeros(len(q))

for i in range(len(q)):
    S[i] = InterpolateSens(log_s, log_q, surv_sens, s[i], q[i], 4, 'cubic')
    #print(S[i])

# setting starting positions for walkers

pos = np.random.rand(50, 4)
nwalkers, ndim = pos.shape

# running emcee

sampler = emcee.EnsembleSampler(nwalkers, ndim, get_log_probability, args = (q_break, sensitivity_data, planet_data, S))
pos_new = sampler.run_mcmc(pos, 400, progress=True)
sampler.reset()
sampler.run_mcmc(pos_new, 6000, progress=True)

flat_samples = sampler.get_chain(flat = True)

# writting results to file
txt = ''
for i in range(ndim):
    results = np.percentile(flat_samples[:, i], [16, 50, 84])
    q = np.diff(results)
    txt += ','+','.join([str(results[1]), str(-q[0]), str(q[1])])

with open("./results/f_params.CSV", 'a', encoding="utf-8") as results_file:
    today = datetime.now()
    results_file.write('\n'+today.strftime("%d/%m/%Y %H:%M:%S")+",Suzuki's" + txt)

# visualization of results

labels = ["A", "m", "n", "p"]

figure = corner.corner(flat_samples, labels=labels,
                       quantiles=[0.16, 0.5, 0.84], bins=40, 
                       show_titles=True, title_kwargs={"fontsize": 12})
plt.savefig("cornerplots/podejscie_suzukiego.png", dpi=600)
