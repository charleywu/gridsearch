#Sample 2D functions from a Gaussian Process RBF prior and save outputs as a json file
#Charley Wu 2017

import numpy as np
import GPy
import json
from matplotlib import pyplot as plt
from sklearn import preprocessing
plt.ioff()

#function to sample values from a gaussian process kernel
def samplePrior(k, imageFileName, xmin=0, xmax=10, normalize = True):
	"""Samples a function from a gaussian process kernel (k), where xmin and xmax specify the square bounds of the 2D function, and normalize=True scales the payoffs between 0 and 1 """
	#Create grid using xmin and xmax
	xx, yy = np.mgrid[xmin:xmax+1, xmin:xmax+1]
	X = np.vstack((xx.flatten(), yy.flatten())).T 
	K = k.K(X) #compute covariance matrix of x X
	s = np.random.multivariate_normal(np.zeros(X.shape[0]), K) #GP prior distribution
	#set min-max range for scaler
	min_max_scaler = preprocessing.MinMaxScaler(feature_range=(0,1))
	#scale to range
	if normalize==True:
		s = min_max_scaler.fit_transform(s)
	#plot figure
	fig = plt.figure()
	plt.matshow(s.reshape(*xx.shape), interpolation='none', cmap=plt.get_cmap('OrRd'))
	plt.xticks(np.arange(0.5,10.5,1), [])
	plt.yticks(np.arange(0.5,10.5,1), [])
	plt.tick_params(axis='both', which='both',length=0)
	plt.grid(color='#0082C1', linestyle='-')
	#save fig
	plt.savefig(imageFileName, bbox_inches='tight', format='pdf')
	plt.close(fig)
	#convert into JSON object
	jsonData = {}
	counter = 0
	for x1 in range(xmin,xmax+1):
		for x2 in range(xmin,xmax+1):
			jsonData[counter] = {'x1':x1, 'x2':x2, 'y':s[counter]}
			counter+=1
	return jsonData, fig


#Create experiment data

#define kernels
kernel1 = GPy.kern.RBF(input_dim=2, variance=1, lengthscale=1)
kernel2 = GPy.kern.RBF(input_dim=2, variance=1, lengthscale=2)
kernelList = [kernel1, kernel2]
filenames = ['experiment2D/kernel1.json', 'experiment2D/kernel2.json']
#number of samples 
samples = 20

#Sample and save output data and plots
for n in range(2):
	outputData = {}
	#sample payout matrix and also create plots
	for sample in range(samples):
		figName = 'experiment2D/images/kernel%i.%i.pdf' % (n, sample)
		(outputData[sample], figure) = samplePrior(kernelList[n], figName)
		
	#save payout matrix
	filename = filenames[n]
	with open(filename, 'w') as fp:
		json.dump(outputData, fp)
