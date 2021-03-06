import numpy as np
from combgam import combgam
import matplotlib.pyplot as plt

x = np.logspace(-6,3,1000)

plind = -1.
log_k = 0.
log_bmin = -6.
alpha = 3.
log_bmax = 1.

ytot, yvals, b = combgam(x, plind, log_k, log_bmin, log_bmax, alpha)
plt.figure(1)
plt.subplot(121)
plt.plot(x, ytot)
plt.yscale('log')
plt.xscale('log')

weights = yvals/np.sum(yvals)

def gamma_sampler(alpha, b, n, plind, log_k):
    sampled = []
    ignored = []
    bs = []
    for i in range(n):
        bstar = np.random.choice(b, p = weights)
        x = np.random.gamma(alpha, bstar)
        #x = np.random.gamma(alpha, 1/bstar)
        if x < 1e-6:
            ignored.append(x)
        else:
            sampled.append(x)
            bs.append(bstar)
    sampled = np.array(sampled)        
    return (sampled, ignored, bs)

xvals, ignoredxvals, bs = gamma_sampler(alpha, b, 5000, plind, log_k)

#print ignoredxvals
plt.subplot(122)
plt.hist(xvals, bins = 10**np.linspace(-6,3,100))
plt.yscale('log')
plt.xscale('log')
plt.show()

#print xvals

datafile_path = "/local/php16lpg/james_code/gamsample1.txt"
datafile_id = open(datafile_path, 'w+')

data = np.array(xvals)

np.savetxt(datafile_id, data, fmt = ['%.10f'])

datafile_id.close()

