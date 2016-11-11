import numpy as np
from combgam import combgam
import numpy.polynomial.polynomial as poly
from statsmodels.distributions.empirical_distribution import ECDF
import matplotlib.pyplot as plt
import scipy.interpolate as interpolate


x = np.logspace(-6,3,1000)
plind = -0.5
log_k = 0.
log_bmin = -3
alpha = 3
ytot, yvals, b = combgam(x, plind, log_k, log_bmin, alpha)
cumy = np.cumsum(ytot/np.sum(ytot))
xvals = []
for i in range(100000):
    u = np.random.uniform(np.min(cumy), np.max(cumy))
    f = interpolate.interp1d(cumy, x)
    xvals.append(f(u))

plt.figure(1)
plt.subplot(121)
plt.plot(x, ytot)
plt.yscale('log')
plt.xscale('log')


plt.subplot(122)
plt.hist(xvals, bins =  10**np.linspace(-6,3,100))
plt.yscale('log')
plt.xscale('log')
plt.show()

