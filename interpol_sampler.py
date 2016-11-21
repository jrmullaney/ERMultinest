import numpy as np
from combgam import combgam
import numpy.polynomial.polynomial as poly
from statsmodels.distributions.empirical_distribution import ECDF
import matplotlib.pyplot as plt
import scipy.interpolate as interpolate


x = np.logspace(-6,3,1000)
plind = 0
log_k = 0.
log_bmin = -6
log_bmax = 1
alpha = 3
ytot, yvals, b = combgam(x, plind, log_k, log_bmin, log_bmax, alpha)
cumy = np.cumsum(ytot/np.sum(ytot))
xvals = []
for i in range(10000):
    u = np.random.uniform(np.min(cumy), np.max(cumy))
    f = interpolate.interp1d(cumy, x)
    xvals.append(f(u))
