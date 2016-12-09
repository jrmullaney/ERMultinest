import numpy as np
import math
import pdb
import matplotlib.pyplot as plt
from qcombgam import qcombgam
import scipy.stats as stats
from norm_function import norm_function

x = np.logspace(-10,5,1000)

log_k = 0.
log_bmin = -10.
log_bmax = 5.
a = 3.
b = np.logspace(log_bmin, log_bmax, 40.)

alpha = 0.06
beta = -1.
log_p = -6.
log_q = 1.

data = np.loadtxt('SSFR_REDD.txt')
ssfr = data[0:1999,0]

#print np.min(ssfr)
#print np.max(ssfr)
#quit()

def gamma_sampler(x, a, log_bmin, log_bmax, n, plind, log_k, log_p, log_q, s):
    ytot, yvals, b = qcombgam(x, plind, log_k, log_bmin, log_bmax, a, log_p, log_q, s)
    weights = yvals/np.sum(yvals)
    sampled = []
    ignored = []
    bs = []
    for i in range(n):
        bstar = np.random.choice(b, p = weights)
        x1 = np.random.gamma(a, bstar)
        #x = np.random.gamma(a, 1/bstar)
        #if x1 < 1e-6:
            #ignored.append(x1)
        #else:
        sampled.append(x1)
        bs.append(bstar)
    sampled = np.array(sampled)        
    return (sampled, ignored, bs)

new_vals = []
plindvals = []
for i in range(ssfr.size):
    plind = alpha * ssfr[i] + beta
    xvals, _, bs = gamma_sampler(x, a, log_bmin, log_bmax, 1, plind, log_k, log_p, log_q, -10.)
    #print xvals
    print i*100/ssfr.size, " % done"
    new_vals.append(xvals)
    plindvals.append(plind)
new_vals = np.hstack(new_vals)
new_vals = np.squeeze(new_vals)

print np.min(plindvals)
print np.max(plindvals)


new_vals = new_vals


#### test ####
#test_plind = alpha * (np.min(ssfr) + (np.max(ssfr) - np.min(ssfr))/2) + beta
#test_vals, _, _ = gamma_sampler(x, a, log_bmin, log_bmax, 10000, test_plind, log_k)

plt.hist(new_vals, bins = 10**np.linspace(-10,3,100))
#plt.hist(test_vals, bins = 10**np.linspace(-10,3,100))
#plt.plot(new_vals, ytot, 'ro')
#plt.plot(ssfr,new_vals, 'x')
plt.xscale('log')
plt.yscale('log')
plt.show()

datafile_path = "/local/php16lpg/james_code/qplind_samp.txt"
datafile_id = open(datafile_path, 'w+')

data = np.array([ssfr, new_vals])
data = data.T

np.savetxt(datafile_id, data, fmt=['%.10f','%.10f'])

datafile_id.close()
