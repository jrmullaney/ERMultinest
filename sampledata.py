import numpy as np
import math
import pdb
import matplotlib.pyplot as plt

#Read in the data
data = np.loadtxt('SSFR_REDD.txt')
ssfr = data[:,0]

#Fix parameter values
alpha = 0.05
beta = 0.1
b = 1
a = alpha*ssfr + beta

#Generate redd from gamma distribution with said parameters

#redd = np.random.gamma(a,b, 10000)

redd = []
for j in a:
    redd1 = np.random.gamma(j,b) 
    redd.append(redd1)
  
#print redd

#Save the txt file

datafile_path = "/local/php16lpg/james_code/SSFR_REDD_samp.txt"
datafile_id = open(datafile_path, 'w+')

data = np.array([ssfr, redd])
data = data.T

np.savetxt(datafile_id, data, fmt=['%.10f','%.10f'])

datafile_id.close()

print np.min(a)
print np.max(a)

plt.hist(np.log(redd), bins = 100)
plt.yscale('log')
plt.show()
