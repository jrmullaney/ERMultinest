import numpy as np
from scipy.stats import gamma
import scipy.special as ss
import matplotlib.pyplot as plt
from norm_function import norm_function
#def dgamma(x, a, b):
#    l = (b**a)/gamma(a)
#    m = x**(a-1.)
#    n = np.exp(-b*x)
#    p = l * m * n

#    o = np.argwhere(np.isnan(p))

#    p[o] = 0

#    return p
#s = 5

def qcombgam(x, plind, log_k, log_bmin, log_bmax, alpha, log_p, log_q, s):

    ytot = 0.*x
    f = np.linspace(0, 19, 20)
    b = 10**(log_bmin + (log_bmax - log_bmin)* 0.05 *f)
    val = norm_function(np.log10(b), log_p, log_q, s)
    #oo = np.where(val <= 1e-05)
    #val[oo] = 0
    #val1 = (np.log10(b) >= log_p)
    #val2 = (np.log10(b) <= log_q)
    #val = val1 * val2
    norm = b**(plind + 1.) * val
    norm = norm / np.sum(norm)
    k = 10**log_k
    yvals = []    

    #plt.figure(1)
    for i in range(b.size):
        #print norm[i]
        y = k*norm[i]*gamma.pdf(x, alpha, 0, b[i])
        o = np.where(y == 0)
        y[o] = 1e-303
        yvals.append(np.max(y))
        ytot = ytot + y
        #plt.plot(x, y)
        #plt.yscale('log')
        #plt.xscale('log')
        #plt.ylim(1e-5, 1e0)
    #plt.plot(x, ytot)
    return (ytot, yvals, b)

"""

x = np.logspace(-12,5,2000)
plind = 1
log_k = 0.
log_bmin = -8.
a = 3.
log_bmax = 3.
log_p = -6.
log_q = 1.
s = 2.

ytot, yvals, b = qcombgam(x, plind, log_k, log_bmin, log_bmax, a, log_p, log_q, s)
#y = x**plind
#plt.plot(x,y)
plt.show()    
  
"""

