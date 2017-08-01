#Sample 1D functions from a Gaussian Process RBF prior and save outputs as a json file
#Charley Wu 2017

import numpy as np
import GPy
import json
from matplotlib import pyplot as plt
from sklearn import preprocessing
plt.ioff()

#function to sample values from a gaussian process kernel
def samplePrior(k, imageFileName, xmin=0, xmax=30, normalize = True):
    """Samples a function from a gaussian process kernel (k), where xmin and xmax specify the square bounds of the 2D function, and normalize=True scales the payoffs between 0 and 1 """
   #Create grid using xmin and xmax
    X = np.arange(xmin,xmax) #X*
    X= X[:,None] #Convert into a matrix
    K = k.K(X,X) #calculate covariance matrix
    s = np.random.multivariate_normal(np.zeros(X.shape[0]), K)
    s = (s-min(s))/(max(s)-min(s))
    fig = plt.figure()
    plt.plot(X,s)
    #save fig
    plt.savefig(imageFileName, bbox_inches='tight', format='pdf')
    plt.close(fig)
    #convert into JSON object
    jsonData = {}
    counter = 0
    for x in range(xmin,xmax):
            jsonData[counter] = {'x':x, 'y':s[counter]}
            counter+=1
    return jsonData, fig


#Create experiment data

#define kernels
kernel1 = GPy.kern.RBF(input_dim=1, variance=1, lengthscale=1)
kernel2 = GPy.kern.RBF(input_dim=1, variance=1, lengthscale=2)
kernelList = [kernel1, kernel2]
filenames = ['experiment1D/kernel1.json', 'experiment1D/kernel2.json']
#number of samples 
samples = 40

#Sample and save output data and plots
for n in range(2):
    outputData = {}
    #sample payout matrix and also create plots
    for sample in range(samples):
        figName = 'experiment1D/images/kernel%i.%i.pdf' % (n, sample)
        (outputData[sample], figure) = samplePrior(kernelList[n], figName)
        
    #save payout matrix
    filename = filenames[n]
    with open(filename, 'w') as fp:
        json.dump(outputData, fp)
