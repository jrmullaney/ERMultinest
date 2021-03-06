import numpy as np
from combgam import combgam
import numpy.polynomial.polynomial as poly
from statsmodels.distributions.empirical_distribution import ECDF
import matplotlib.pyplot as plt
import scipy.interpolate as interpolate

x = np.logspace(-6,3,1000)

plind = 0.
log_k = 0.
log_bmin = -6.
alpha = 3
log_bmax = 1

ytot, yvals, b = combgam(x, plind, log_k, log_bmin, log_max, alpha)



#ytote = np.copy(ytot)
#o = np.where(ytot <= 1e-10)
#ytot[o] = 0

#ysums = ytot/np.sum(ytot)
#ecdf = np.cumsum(ysums)

#plt.plot(x, ecdf)
plt.figure(1)
plt.subplot(131)
plt.plot(x, ytot)
plt.yscale('log')
plt.xscale('log')




#def sampler(n):
#    for i in range(n):
#        u = np.random.uniform()
#        x =


def sampler(x, weights, n):
    sampled = []
    for i in range(n):
        sampled.append(np.random.choice(x, p = weights))
    return sampled
                       
weights = ytot/np.sum(ytot)


vals = sampler(x, weights, 10000)
plt.subplot(132)
plt.hist(vals, bins = 10**np.linspace(-6,3,100))
plt.yscale('log')
plt.xscale('log')


cumy = np.cumsum(ytot/np.sum(ytot))
plt.subplot(133)
plt.plot(x, cumy)
plt.yscale('log')
plt.xscale('log')


plt.show()



