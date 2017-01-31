import pymultinest
#import math
import numpy as np
#import pdb
import numpy.polynomial.polynomial as poly
from scipy.integrate import quad
from scipy.stats import gamma
from scipy.stats import loggamma
import scipy.special as special

log_bmax = 1.
log_bmin = -6.
a = 1.
f = np.linspace(0, 19, 20)
b = np.array(10**(-10. + 15. * 0.05 * f))

def myprior(cube, ndim, nparams):
     
    #k is normalisation
    #beta is PL index
    #bmin is low turnover
    #bmax is high turnover
    log_k = -2 + 4.*cube[0]  #Normalisation
    beta = -4. + 8.*cube[1]  #PL index
    log_bmin = -10 + 7*cube[2] #low turnover
    log_bmax = -3 + 6*cube[3] #low turnover
          
    cube[0] = log_k 
    cube[1] = beta
    cube[2] = log_bmin
    cube[3] = log_bmax    
    
def myloglike(cube, ndim, nparams):

    #Model Parameters
    k = 10.**cube[0]
    beta = cube[1]
    log_bmin = cube[2]
    log_bmax = cube[3]
    
    #Calculating the norm pdfs/cdfs for each galaxy
    pnorm2d = pb2d ** beta
    cnorm2d = cb2d ** beta
    pnorm2d[np.log10(pb2d) > log_bmax] = 0.
    pnorm2d[np.log10(pb2d) < log_bmin] = 0.
    
    #Ensure the normalisations sum to 1:
    pn = np.sum(pnorm2d, axis=0)
    cn = np.sum(cnorm2d, axis=0)
    pn = np.tile(pn[np.newaxis,:],(b.size,1))
    cn = np.tile(cn[np.newaxis,:],(b.size,1))
    pnorm2d = pnorm2d/pn
    cnorm2d = cnorm2d/cn
    
    #Multiply the pdfs and cdfs by the norms and sum
    #to get total likelihoods.
    norm_pdf = pnorm2d * pdf2d
    sum_norm_pdf = np.sum(norm_pdf, axis = 0)        
    
    norm_cdf =  cnorm2d * cdf2d
    sum_norm_cdf = np.sum(norm_cdf, axis = 0)
    
    #First part of log lik
    lik1 = np.log(k * sum_norm_pdf)
    
    #For the second part, it's the sum of the integrals:
    lik2 = k * np.sum(sum_norm_cdf)
    
    #Likelihood
    ln_l = np.sum(lik1) - lik2
    return ln_l

#Name and number of parameters our problem has:
parameters = ['k', 'beta', 'bmin', 'bmax']
n_params = len(parameters)

#Read the sampled data
data = np.loadtxt('plind_samp.txt')
data = data[0:2000,:]
sam = data[:,1]
nsam = sam.size

#Detected sample:
UL = 1e-5*np.ones(nsam)
o = np.where(sam >= UL)
redd = sam[o]
ndet = redd.size

#This calculates the lik of each detected redd for each gamma:
pdf2d = np.zeros((b.size, ndet))
for i in range(ndet):
    pdf_val = gamma.pdf(redd[i], a, 0, b)
    pdf2d[:,i] = pdf_val

#This calculates the CDF for all the ULs for each gamma:
cdf2d = np.zeros((b.size, nsam))
for i in range(nsam):
    cdf_val = 1. - gamma.cdf(UL[i], a, 0, b)
    cdf2d[:,i] = cdf_val

pb2d = np.tile(b[:,np.newaxis], (1, ndet))
cb2d = np.tile(b[:,np.newaxis], (1, nsam))

#ssfr2d = np.tile(ssfr[np.newaxis,:], (b.size, 1))
#sam_ssfr2d = np.tile(sam_ssfr[np.newaxis,:], (b.size, 1))

#Run MultiNest
pymultinest.run(myloglike, myprior, n_params, importance_nested_sampling = False, resume = False, verbose = True, sampling_efficiency = 'model', n_live_points = 500, outputfiles_basename='chains/PL_twogal_redd-')
#Analyse the results and plot the marginals
ana = pymultinest.Analyzer(n_params = n_params, outputfiles_basename='chains/PL_twogal_redd-')

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

plt.savefig("chains/marginals_multinest_twogal_redd.pdf") #, bbox_inches='tight')
#show("chains/marginals_multinest.pdf")
