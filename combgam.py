import numpy as np
from scipy.stats import gamma
import scipy.special as ss
import matplotlib.pyplot as plt

#def dgamma(x, a, b):
#    l = (b**a)/gamma(a)
#    m = x**(a-1.)
#    n = np.exp(-b*x)
#    p = l * m * n

#    o = np.argwhere(np.isnan(p))

#    p[o] = 0

#    return p


def combgam(x, plind, log_k, log_bmin, log_bmax, alpha):

    ytot = 0.*x
    f = np.linspace(0, 19, 20)
    b = 10**(log_bmin + (log_bmax - log_bmin)* 0.05 *f)
    norm = b**(plind +1.)
    k = 10**log_k
    yvals = []
    

    #plt.figure(1)
    for i in range(b.size):
        #print norm[i]
        y = k*norm[i]*gamma.pdf(x, alpha, 0, b[i])
        o = np.where(y <= 1e-10)
        y[o] = 0
        yvals.append(np.max(y))
        ytot = ytot + y
        #plt.plot(x, y)
        #plt.yscale('log')
        #plt.xscale('log')
        #plt.ylim(1e-5, 1e3)
    #plt.plot(x, ytot)
    return (ytot, yvals, b) 

#plind = 0.
#log_k = 0.
#log_bmin = -6.
#alpha = 3

#ytot, yvals, b = combgam(x, plind, log_k, log_bmin, alpha)
#print yvals
#y = x**plind
#plt.plot(x,y)
#plt.show()    
    
