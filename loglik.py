import numpy as np
from scipy.special import gamma
from scipy.special import gammainc
import matplotlib.pyplot as plt

def genb(alpha, bmin):
    bmax = np.log10(bmin) + 0.3*19.
    
    b = np.logspace(np.log10(bmin), bmax, 20)

    return b

def cgamma(x, bmin, alpha, plind):
    
    b = genb(alpha, bmin)
    norm = b**(-plind)
    bt = b[:,np.newaxis]
    normt = norm[:,np.newaxis]

    l = (bt**alpha)/gamma(alpha)
    m = x**alpha
    n = np.exp(-bt*x)
    
    sep = normt*l*m*n

    tot = np.sum(sep, axis=0)

    #for curves in sep[0,:]
    #    plt.plot(x, curves)
    #plt.show()
    
    return tot

def intgamma(xmin, bmin, alpha, plind):
    
    b = genb(alpha, bmin)
    norm = b**(-plind)

    #Scipy's igamma is (1/gamma(a))*int(x^...)
    integ = norm*(1.-gammainc(alpha,b*xmin))

    #Sum the integrand from all of the separate gamma
    #distributions:
    return np.sum(integ) 
