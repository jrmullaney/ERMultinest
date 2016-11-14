import pymultinest
import math
import numpy as np
import pdb
import numpy.polynomial.polynomial as poly
from scipy.integrate import quad
from scipy.stats import gamma
from scipy.stats import loggamma
import scipy.special as special

log_bmax = 1.

def myprior(cube, ndim, nparams):
    plind = 2.*cube[0] - 1. 
    log_k = 12.*cube[1] - 6.
    log_bmin = -12.*cube[2] 
    a = 6.*cube[3] 
    
    k = 10.**log_k
    f = np.linspace(0, 19, 20)
    b = 10**(log_bmin + (log_bmax - log_bmin) * 0.05 * f)
    norm = b**(plind)

    
    cube[0] = plind
    cube[1] = k
    cube[2] = log_bmin
    cube[3] = a

def myloglike(cube, ndim, nparams):

    #Model Parameters
    plind = cube[0]
    k = cube[1]
    log_bmin = cube[2]
    a = cube[3]
    
    f = np.linspace(0, 19, 20)
    nsam = redd.size
    b = 10.**(log_bmin + (log_bmax - log_bmin) * 0.05 * f)     
    norm = b**(plind)
    UL = 1e-3
    
    #Detected sample:
    o = np.where(redd >= UL)
    det = redd[o]
    ndet = det.size
    lndet = np.log(det)
    av_lndet = np.mean(lndet)

    #First part of likelihood function:
    lik = []
    pdf = []
    pdf1 = np.zeros( (det.size, b.size) )
    #calculate the pdf of 1 det, for all bs
    
    for j in range(b.size):
        pdf1[:,j] = k * norm[j] * gamma.pdf(det, a, 0, b[j])
        
    lik1 = np.sum(pdf1, axis = 1)
    loglik1 = np.log(lik1)
    loglikT = np.sum(loglik1)
    
            #a = ndet*np.log(k) + np.log(np.sum(norm[i])) + np.log(np.sum(gamma.pdf(det, a, 0, b[i]))))
        #lik.append(a)
           
    
      
    
    #For the second part, need to integrate:
    cdf = []
    #cdf_est = 1 - (ndet / nsam)
    for i in range(b.size):
        #cdf = 1 - (ligf(a, b[i]*1e-3)/math.gamma(a))
        cdfval = (norm[i]*(1 - gamma.cdf(UL, a, 0, b[i])))
        cdf.append(cdfval)
    
    int2 = np.sum(cdf)
    lik2 = k*int2
    
    
    #Likelihood
    ln_l = np.sum(loglikT - (nsam*lik2))
    return ln_l

#Name and number of parameters our problem has:
parameters = ['plind', 'log_k', 'log_bmin', 'a']
n_params = len(parameters)

#Read the sampled data
data = np.loadtxt('gamsample1.txt')
redd = data

#Run MultiNest
pymultinest.run(myloglike, myprior, n_params, importance_nested_sampling = False, resume = False, verbose = True, sampling_efficiency = 'model', n_live_points = 500, outputfiles_basename='chains/PL_redd-')

#Analyse the results and plot the marginals
ana = pymultinest.Analyzer(n_params = n_params, outputfiles_basename='chains/PL_redd-')

import matplotlib.pyplot as plt
plt.clf()

# Here we will plot all the marginals and whatnot, just to show off
# You may configure the format of the output here, or in matplotlibrc
# All pymultinest does is filling in the data of the plot.

# Copy and edit this file, and play with it.

p = pymultinest.PlotMarginalModes(ana)
plt.figure(figsize=(5*n_params, 5*n_params))
#plt.subplots_adjust(wspace=0, hspace=0)
for i in range(n_params):
	plt.subplot(n_params, n_params, n_params * i + i + 1)
	p.plot_marginal(i, with_ellipses = True, with_points = False, grid_points=50)
	plt.ylabel("Probability")
	plt.xlabel(parameters[i])
	
	for j in range(i):
		plt.subplot(n_params, n_params, n_params * j + i + 1)
		#plt.subplots_adjust(left=0, bottom=0, right=0, top=0, wspace=0, hspace=0)
		p.plot_conditional(i, j, with_ellipses = False, with_points = True, grid_points=30)
		plt.xlabel(parameters[i])
		plt.ylabel(parameters[j])

plt.savefig("chains/marginals_multinest_redd.pdf") #, bbox_inches='tight')
#show("chains/marginals_multinest.pdf")
