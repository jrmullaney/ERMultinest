import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt



#x = np.linspace(-10, 10, 1000)


def norm_function(x, x1, x2, s):
    pdfvals = []
    for i in x:
        if (i < x1):
            pdf = 1.02 / (1 + np.exp(-s*(i-(x1 - 2))))
            #pdf = 0 
            pdfvals.append(pdf)
        elif (x1 <= i) and (i < x2):
            pdf = 1
            pdfvals.append(pdf)
        else:
            pdf =  1.02 / (1 + np.exp(s*(i-(x2 + 2))))
            #pdf = 0 
            pdfvals.append(pdf)
    pdfvals = np.array(pdfvals)
    o = np.where(pdfvals > 1)
    pdfvals[o] = 1
    return pdfvals

"""

#log_x1 = -6
#log_x2 = 1
#x1 = 10. ** log_x1
#x2 = 10. ** log_x2
s = 100
x1 = -5
x2 = 2
vals = norm_function(x, x1, x2, s)


plt.plot(x, vals)
#plt.xscale('log')
#plt.yscale('log')
plt.ylim(0,1.2)
plt.show()

"""
