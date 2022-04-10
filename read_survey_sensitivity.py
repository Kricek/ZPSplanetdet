"""
Read survey sensitivity data and plot the chart
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

sensitivity_data = pd.read_csv('data/survey_sensitivity2.dat', header=None)
data = sensitivity_data.to_numpy()
log_s = data[0,2:]
log_q = data[1:,1]
surv_sens = data[1:,2:]

plt.yscale("log")
plt.xscale("log")
cs = plt.contour(10**log_s, 10**log_q, surv_sens, [5,20,50,200,500,700], colors = 'k')
plt.clabel(cs)
ax = plt.gca()
ax.set_title("Survey Sensitivity, $S(\\log(s), \\log(q))$")
ax.set_xlabel("Separation, $s$ $(\\theta_E)$")
ax.set_ylabel("Mass ratio, $q$")
ax.set_xticks([0.2,0.5,1,2,5])
ax.set_xticklabels([0.2,0.5,1,2,5])
plt.savefig("survey_sensitivity.png")
