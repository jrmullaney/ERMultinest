import numpy as np
import math
import pdb
import matplotlib.pyplot as plt
from combgam import combgam
from gamma_sampler import gamma_sampler


x = np.logspace(-6,3,1000)
plind = 0.
log_k = 0.
log_bmin = -6.
a = 3
log_bmax = 1
alpha = 1
beta = -0.25
data = np.loadtxt('SSFR_REDD.txt')
ssfr = data[:,0]


ytot, yvals, b = combgam(x, plind, log_k, log_bmin, log_bmax, alpha)
xvals, ignoredxvals, bs = gamma_sampler(a, b, 5000, plind, log_k)

bs = np.array(bs)
xvals = np.array(xvals)
ssfr = ssfr[0:xvals.size]


new_vals = []
for i in range(bs.size):
    plind = alpha * ssfr[i] + beta
    newval = bs[i] ** plind * xvals[i]
    new_vals.append(newval)


datafile_path = "/local/php16lpg/james_code/REDD_SSFR_plind_samp.txt"
datafile_id = open(datafile_path, 'w+')

data = np.array([ssfr, new_vals])
data = data.T

np.savetxt(datafile_id, data, fmt = ['%.10f', '%.10f'])

datafile_id.close()

