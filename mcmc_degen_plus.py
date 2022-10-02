from cmath import nan
import numpy as np
import emcee
import pandas as pd
import corner
import matplotlib.pyplot as plt
from tqdm import tqdm
from datetime import date

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
    A, m, n, p, x_4, x_19, x_22, x_24 = params
    model_params = (A, m, n, p)
    i = 0
    
    log_L0 = -1 * get_N_exp(model_params, q_break, sensitivity_data)
    
    while i < len(q):        
        if(w[i] == 1):    
            log_L0 += np.log(A) + np.log(get_massratio(model_params, q_break, s[i], q[i])) + np.log(S[i])
        i += 1
       
    if(x_4 < w[4]):
        log_L0 += np.log(A) + np.log(get_massratio(model_params, q_break, s[4], q[4])) + np.log(S[4])
    else:
        log_L0 += np.log(A) + np.log(get_massratio(model_params, q_break, s[5], q[5])) + np.log(S[5])

    if(x_19 < w[19]):
        log_L0 += np.log(A) + np.log(get_massratio(model_params, q_break, s[19], q[19])) + np.log(S[19])
    else:
        log_L0 += np.log(A) + np.log(get_massratio(model_params, q_break, s[20], q[20])) + np.log(S[20])    

    if(x_22 < w[22]):
        log_L0 += np.log(A) + np.log(get_massratio(model_params, q_break, s[22], q[22])) + np.log(S[22])
    else:
        log_L0 += np.log(A) + np.log(get_massratio(model_params, q_break, s[23], q[23])) + np.log(S[23])

    if(x_24 < w[24]):
        log_L0 += np.log(A) + np.log(get_massratio(model_params, q_break, s[24], q[24])) + np.log(S[24])
    else:
        log_L0 += np.log(A) + np.log(get_massratio(model_params, q_break, s[25], q[25])) + np.log(S[25])

    return log_L0

def get_log_prior(params):
    """
    uninformative prior pdf
    """
    A, m, n, p, x_4, x_19, x_22, x_24 = params 
    if  A > 0 and x_4 > 0 and x_4 < 1 and x_19 > 0 and x_19 < 1 and x_22 > 0 and x_22 < 1 and x_24 > 0 and x_24 < 1:
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

x_4 = 0.5
x_19 = 0.5
x_22 = 0.5
x_24 = 0.5

params = [A, m, n, p, x_4, x_19, x_22, x_24]

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

pos = np.random.rand(64, 8)
nwalkers, ndim = pos.shape

# running emcee

sampler = emcee.EnsembleSampler(nwalkers, ndim, get_log_probability, args = (q_break, sensitivity_data, planet_data, S))
pos_new = sampler.run_mcmc(pos, 400, progress=True)
sampler.reset()
sampler.run_mcmc(pos_new, 6000, progress=True)


# obtaining the flat samples

flat_samples = sampler.get_chain(flat = True)
A_plot, m_plot, n_plot, p_plot, x4_plot, x19_plot, x22_plot, x24_plot = flat_samples.T

# writting results to file
txt = ''
for i in range(ndim):
    results = np.percentile(flat_samples[:, i], [16, 50, 84])
    q = np.diff(results)
    txt += ','+','.join([str(results[1]), str(-q[0]), str(q[1])])

with open("./results/f_params.CSV", 'a', encoding="utf-8") as results_file:
    today = date.today()
    results_file.write('\n'+today.strftime("%d/%m/%Y %H:%M:%S")+",new (no deviation)" + txt)

# printing the fraction of values above/below the limit

filtr1 = (x4_plot < w[5])
filtr2 = (x19_plot < w[20])
filtr3 = (x22_plot < w[22])
filtr4 = (x24_plot < w[24])
filtry = [filtr1, filtr2, filtr3, filtr4]

x4_below = x4_plot[filtr1]
x4_above = x4_plot[~filtr1]

x19_below = x19_plot[filtr2]
x19_above = x19_plot[~filtr2]

x22_below = x22_plot[filtr3]
x22_above = x22_plot[~filtr3]

x24_below = x24_plot[filtr4]
x24_above = x24_plot[~filtr4]

print("n(x4<w4)/n: ", np.sum(filtr1)/len(filtr1))

print("n(x19<w19)/n: ", np.sum(filtr2)/len(filtr2))

print("n(x22<w22)/n: ", np.sum(filtr3)/len(filtr3))

print("n(x24<w4)/n: ", np.sum(filtr4)/len(filtr4))

# visualization of the results

labels1 = ["A", "m", "n", "p"]
labels2 = ["x4_below"]
labels3 = ["x4_above"]

results = np.array([A_plot, m_plot, n_plot, p_plot])#_below, x4_above])
results = results.T

figure = corner.corner(results, labels=labels1, quantiles=[0.16, 0.5, 0.84], bins=40, show_titles=True, title_kwargs={"fontsize": 12})
plt.savefig("cornerplots/nowe_podejscie.png", dpi=600)
plt.clf()

names = ["x_4", "x_19", "x_22", "x_24"]

for i in range(4):

    sum1 = np.sum(filtry[i])
    sum2 = np.sum(~filtry[i])

    results_below = np.array([A_plot[filtry[i]], m_plot[filtry[i]], n_plot[filtry[i]], p_plot[filtry[i]]])
    results_above = np.array([A_plot[~filtry[i]], m_plot[~filtry[i]], n_plot[~filtry[i]], p_plot[~filtry[i]]])

    figure = corner.corner(results_below.T, quantiles=[0.16, 0.5, 0.84], labels=labels1, bins=40, color="blue", label = "x4_below")                                   
    figure = corner.corner(results_above.T, quantiles=[0.16, 0.5, 0.84], labels=labels1, bins=40, color="red", label = "x4_above", weights=(sum1/sum2)*np.ones(sum2), fig=figure)

    plt.savefig("cornerplots/" + names[i] + "/" + names[i] + "_above_below.png", dpi=600)                       
    plt.clf()

    figure = corner.corner(results_below.T, quantiles=[0.16, 0.5, 0.84], show_titles=True, labels=labels1, color="blue", label = "x4_below")
    plt.savefig("cornerplots/" + names[i] + "/" + names[i] + "_below.png")                          
    figure = corner.corner(results_above.T, quantiles=[0.16, 0.5, 0.84], show_titles=True, labels=labels1, color="red", label = "x4_above")
    plt.savefig("cornerplots/" + names[i] + "/" + names[i] + "_above.png")


