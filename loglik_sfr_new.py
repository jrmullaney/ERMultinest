import pymultinest
import math
import numpy as np
import pdb
import numpy.polynomial.polynomial as poly
from scipy.integrate import quad

def myprior(cube, ndim, nparams):
    alpha = 0.1*cube[0]
    log_k = 10.*cube[1] - 5
    beta = 0.2*cube[2] 
    #b = 10.*cube[3]-5.
    #b = 1
    k = 10.**log_k
    
    cube[0] = alpha
    cube[1] = k
    cube[2] = beta
    #cube[3] = b
    

def integrand(x, a):
    
    return (x**(a-1.))*np.exp(-x)
    
def myloglike(cube, ndim, nparams):

    #Model parameters:
    alpha = cube[0]
    k = cube[1]
    beta = cube[2]
    #b = cube[3]
    #b = 1
    #Characteristics of data:
    #Full sample:
    nsam = redd.size

    #Detected sample:
    o = np.where(redd >= 1e-3)
    det = redd[o]
    ssfrdet = ssfr[o]
    ndet = det.size
    lndet = np.log(det)
    av_lndet = np.mean(lndet)

    #First part of likelihood function:
    #lik1 = np.sum(np.log(k) + alpha*lndet - det)

    a = beta + alpha*ssfrdet

    #alndet = [x*y for x, y in zip(a, lndet)]

    lik1 = ndet*np.log(k) + np.sum(a*lndet) - np.sum(det)
    
    #For the second part, need to integrate:
    #int2 = abs(quad(integrand, 1e-3, np.inf, args = a)[0])
    int2 = np.sum(10.**poly.polyval(a, pars))
    #print(int/int2)
    lik2 = k*int2
    #Log-likelihood
    ln_l = lik1 - nsam*lik2
    return ln_l

#Name and number of parameters our problem has:
parameters = ["alpha", 'log_k', 'beta']
n_params = len(parameters)

#Read in the data:
data = np.loadtxt('SSFR_REDD_samp.txt')
ssfr = data[:,0]
redd = data[:,1]
nsam = redd.size

#Create a polynomial fit for a:
avals = np.linspace(-4,4,100)
int = np.zeros(avals.size)
ind = 0
for i in avals:
    qval = quad(integrand, 1e-3, np.inf, args=i)
    int[ind] = np.log10(abs(qval[0]))
    ind += 1
pars = poly.polyfit(avals, int, 12)

#Run MultiNest
pymultinest.run(myloglike, myprior, n_params, importance_nested_sampling = False, resume = False, verbose = True, sampling_efficiency = 'model', n_live_points = 500, outputfiles_basename='chains/PL_SFR_new-')

#Analyse the results and plot the marginals
ana = pymultinest.Analyzer(n_params = n_params, outputfiles_basename='chains/PL_SFR_new-')

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

plt.savefig("chains/marginals_multinest_new.pdf") #, bbox_inches='tight')
#show("chains/marginals_multinest.pdf")

