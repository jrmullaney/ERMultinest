import pymultinest
import math
import numpy as np
import pdb
import numpy.polynomial.polynomial as poly
from scipy.integrate import quad
from scipy.stats import gamma
from scipy.stats import loggamma
import scipy.special as special
from norm_function import norm_function


log_bmax = 2.
log_bmin = -10.
a = 1.
f = np.linspace(0, 19, 20)
b = np.array(10**(log_bmin + (log_bmax - log_bmin) * 0.05 * f))
s = 5

def myprior(cube, ndim, nparams):
     
    #I'm saying that p is bmin, q is bmax:
    log_k = 10 * cube[0] - 5.
    alpha = 4. * cube[1] - 2.  
    beta = 4. * cube[2] - 2.
    log_p = 3.*cube[3] - 7.
    log_q = 4. * cube[4] - 1.
    #s = cube[5]*10
    
    cube[0] = log_k
    cube[1] = alpha
    cube[2] = beta
    cube[3] = log_p
    cube[4] = log_q
    #cube[5] = s
    
def myloglike(cube, ndim, nparams):

    #Model Parameters
    
    log_k = cube[0]
    alpha = cube[1]
    beta = cube[2]
    log_p = cube[3]
    log_q = cube[4] 
    #s = cube[5]
    
    k = 10.** log_k
    p = 10.** log_p
    q = 10.** log_q
    
    #Subsetting according to parameters
    pval2d = np.tile(val[:,np.newaxis], (1, ndet))
    cval2d = np.tile(val[:,np.newaxis], (1, nsam))
   
    #Calculating the norm pdfs/cdfs for each galaxy
    pplind = alpha * ssfr2d + beta
    cplind = alpha * sam_ssfr2d + beta
    
    #Normalise to give PL index:
    pnorm2d = (pb2d ** (pplind))
    cnorm2d = (cb2d ** (cplind))

    #Normalise to impose bmin and bmax:
    pmm2d = (1./(1.+np.exp(-10.*(pb2d-p))))*\
            (1./(1.+np.exp(-k*(pb2d-q))))
    cmm2d = (1./(1.+np.exp(-10.*(cb2d-p))))*\
            (1./(1.+np.exp(-k*(cb2d-q))))
    pnorm2d = pmm2d * pnorm2d
    cnorm2d = cmm2d * cnorm2d
    
    #Normalise to 1:
    pnorm2d = pnorm2d / np.sum(pnorm2d, axis = 0)
    cnorm2d = cnorm2d / np.sum(cnorm2d, axis = 0)

    #Apply normalisations to pdfs:
    norm_pdf = pnorm2d * pdf2d
    norm_cdf = cnorm2d * cdf2d
    
    #Sum over the b's to give probability for each Edd:
    sum_norm_pdf = np.sum(norm_pdf, axis = 0)
    sum_norm_cdf = np.sum(norm_cdf, axis = 0)
    
    #First part of log lik
    lik1 = np.log(k * sum_norm_pdf)
    
    #econd part of log lik:
    lik2 = k * np.sum(sum_norm_cdf)
          
    #Likelihood
    ln_l = np.sum(lik1) - lik2

    return ln_l

#Name and number of parameters our problem has:
parameters = ['k', 'alpha', 'beta', 'log_p', 'log_q']
n_params = len(parameters)

#Read the sampled data
data = np.loadtxt('qplind_samp.txt')
data = data[0:500,:]
sam_ssfr = data[:,0]
nsam = sam_ssfr.size
sam_redd = data[:,1]

#Detected sample:
UL = 1e-10
o = np.where(data[:,1] >= UL)
det = data[o]
redd = det[:,1]
ssfr = det[:,0]
ndet = ssfr.size

pdf2d = np.zeros((b.size, redd.size))
for i in range(redd.size):
    pdf_val = gamma.pdf(redd[i], a, 0, b)
    pdf2d[:,i] = pdf_val

#pdfsum = []
#for i in range(redd.size):
#    pdfval = np.sum(gamma.pdf(redd[i], a, 0, b))
#    pdfsum.append(pdfval)

cdf = 1 - gamma.cdf(UL, a, 0, b)
cdf2d = np.tile(cdf[:,np.newaxis], (1, nsam))

pb2d = np.tile(b[:,np.newaxis], (1, ndet))

cb2d = np.tile(b[:,np.newaxis], (1, nsam))

ssfr2d = np.tile(ssfr[np.newaxis,:], (b.size, 1))

sam_ssfr2d = np.tile(sam_ssfr[np.newaxis,:], (b.size, 1))



#Run MultiNest
pymultinest.run(myloglike, myprior, n_params, importance_nested_sampling = False, resume = False, verbose = True, sampling_efficiency = 'model', n_live_points = 500, outputfiles_basename='chains/q_PL_twogal_redd_ssfr-')
#Analyse the results and plot the marginals
ana = pymultinest.Analyzer(n_params = n_params, outputfiles_basename='chains/q_PL_twogal_redd_ssfr-')

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

plt.savefig("chains/q_marginals_multinest_twogal_redd_ssfr.pdf") #, bbox_inches='tight')
#show("chains/marginals_multinest.pdf")
