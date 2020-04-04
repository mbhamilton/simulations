'''
Created on September 20, 2012
Code ported from Matlab code authored by Dr Matthew Hamilton of Georgetown University
The code below draws the H-W parabola on a definetti diagram and plots points that represent
a set of genotype frequencies for one locus with two alleles
@author: marieidowu
'''
from __future__ import division #using python 3.0 behaviour of division which will always return a float instead of casting to an int
import matplotlib.figure as fig
import matplotlib.lines as ln
import numpy as np
import random as rand
from math import sqrt, floor
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import logging as lg

def plot(model,clear,time, mu=None, nu=None):
	
	#create a Figure instance 
	figure = fig.Figure((12,10))
	figure.set_facecolor('#ececec')
	#use the Figure to create one or more Axes or Subplot instances
	#axes = figure.add_axes([1, time, -0.01,1.01])
	axes = figure.add_subplot(1,1,1)
	
	#set axes labels 
	figure.suptitle('Mutation Simulations', size=18,color='#003E7E',family='serif', weight='bold')
	
	if not clear:
		print 'not clear'
		pizero = [0.05,0.25,0.5,0.75,0.95]
		numberOfFrequencies = np.size(pizero)
		frequenciesOverTime = np.zeros((numberOfFrequencies, time+1))
		frequenciesOverTime[:,0] = pizero
		
		
		if model == 0:
			plot_irreversible_mutation(time, mu, frequenciesOverTime)
	
		elif model ==1:
			plot_twoway_mutaion(time, mu, nu,frequenciesOverTime)
		
		format = ['b-','r-','g-','m-','c-']
		for i in xrange(numberOfFrequencies):
			axes.plot(frequenciesOverTime[i,:],format[i],label=r'$p_{0}\ =$'+str(pizero[i]))
			
		axes.set_ylabel("Frequency of A allele (p)", fontdict={'size':16,'color':'#002147','family':'serif'})
		axes.set_xlabel("Generations", fontdict={'size':16,'color':'#002147','family':'serif'})
		axes.set_xlim((1,time))
		axes.set_ylim((-0.01,1.01))
		 #set a legend for the plot using the labels specified in the above
		#plot commands, set it in the best location with fontsize of 12
		axes.legend(loc='best',prop={'size':14})
	#place figure on Agg FigureCanvas
	canvas = FigureCanvas(figure)

	return canvas


def plot_irreversible_mutation(time, mu,frequenciesOverTime):
	for i in xrange(1,time+1):
		frequenciesOverTime[:,i] = ((1-mu)**i)*frequenciesOverTime[:,0] #* is element-wise multiplication, dot() matrix multiplication for arrays 

def plot_twoway_mutaion(time, mu, nu, frequenciesOverTime):
	"""
	calculate the allele frequency with one-way mutation
	mutation is reversible: A -> a and a -> A 
	"""
	for i in xrange(1,time+1):
		frequenciesOverTime[:,i] = (1-mu)*frequenciesOverTime[:,i-1] + (1-frequenciesOverTime[:,i-1])*nu #* is element-wise multiplication, dot() matrix multiplication for arrays 
		
def plot_neutral_mutation(Ne, gens, interval, clear):
	"""
	calculate the allele frequency with one-way mutation mutation is from A allele (frequency p) to the a allele (frequency 1-p)
  
	 graph the frequency of alleles over time that that initially appeared as new mutations with a frequency of 1/(2Ne)
	 """
	 #create a Figure instance 
	figure = fig.Figure((12,10))
	figure.set_facecolor('#ececec')
	#use the Figure to create one or more Axes or Subplot instances
	#axes = figure.add_axes([1, time, -0.01,1.01])
	axes = figure.add_subplot(1,1,1)
	#set axes labels 
	figure.suptitle('Mutation Simulations', size=18,color='#003E7E',family='serif', weight='bold')
	axes.set_ylabel("Allele frequency", fontdict={'size':16,'color':'#002147','family':'serif'})
	axes.set_xlabel("Generation", fontdict={'size':16,'color':'#002147','family':'serif'})
	axes.set_xlim((1,gens))
	axes.set_ylim((0.0,1.05))
	
	if not clear:
		
		numMuts = int(floor(gens/interval))
		frequencies = np.zeros((numMuts,gens))
		mutationCount = 0 #(1 in matlab) counter for mutations introduced
		mutation = 1 / (2 * Ne)
		frequencies[0,0] = mutation #introduce first mutation
		for i in xrange(1,gens): #range is up to but not including end number i.e if gens=5 [2,3,4]
			frequencies = drift(frequencies, mutationCount, Ne,i) #allele frequencies experience genetic drift
		#	 print 'frequencies matrix ' ,frequencies
			#print 'mod of %d and interval %d is %d' %(i,interval,(i%interval))
			if((i%interval) == 0): #another genertic drift
				if(mutationCount < numMuts-1):
					mutationCount +=1
				frequencies[mutationCount,i] =mutation
			sumOfFrequencies = np.sum(frequencies,0)
		ancestralFrequency = 1-sumOfFrequencies
		
		colours = ['k','m','r','g','b']
		lengthOfColours = len(colours)
		
		for m in xrange(numMuts): #use numMuts rather than mutationCount numMut=num of rows in frequencies
			index = m%lengthOfColours
			colsym =colours[index] 
			axes.plot(frequencies[m,:],colsym)
		
	canvas = FigureCanvas(figure)

	return canvas

def drift(pmatrix, nalleles, N, gen):
	'''
	allele freqs drift, uses Kimura's psuedo random variable method.
	pmatrix is the original matrix of allele frequencies (columns)and for each allele (rows). 
	'''
	#print 'pmatrix: %s' %p_matrix
	drift_matrix = pmatrix.copy()
	
	for k in xrange(nalleles+1): #count from 0 not 1
		if (pmatrix[k,gen-1] < 1.0) & (pmatrix[k,gen-1] > 0.0): # may need to be gen -2
			delta = (2*rand.random()-1)*sqrt((3*pmatrix[k,gen-1]*(1-pmatrix[k,gen-1]))/(2*N))
			if pmatrix[k,gen-1]+delta > 1.0:
				drift_matrix[k,gen] =1.0
			elif pmatrix[k,gen-1]+delta < 0.0:
				drift_matrix[k,gen] =0.0
			else:
				drift_matrix[k,gen] = pmatrix[k,gen-1] +delta
		else:
			drift_matrix[k,gen] =pmatrix[k,gen-1] #allele is no longer segregating so freq is constant through time
	#print 'drift_matrix: %s' %drift_matrix
	return drift_matrix	 

def validate_mu(mu):
	" Checks that the value of mu is between 0.01 and 0.00000001"
	min_mu = 0.00000001
	max_mu = 0.01
	if mu < min_mu:
		mu = min_mu
	elif mu > max_mu:
		mu = max_mu
	return mu

def validate_ne(ne):
	" Checks that the value of Ne is between 5 and 10000"
	min_ne = 5
	max_ne = 10000
	if ne < min_ne:
		ne = min_ne
	elif ne > max_ne:
		ne = max_ne
	return ne
		
def validate_generations_between(gensbetween):
	" Checks that the value of Ne is between 5 and 10000"
	min_gensbetween = 5
	max_gensbetween = 50
	if gensbetween < min_gensbetween:
		gensbetween = min_gensbetween
	elif gensbetween > max_gensbetween:
		gensbetween = max_gensbetween
	return gensbetween

def validate_total_generations(totalgens):
	" Checks that the value of total generations is between 100 and 1000"
	min_totalgens = 100
	max_totalgens = 1000
	if totalgens < min_totalgens:
		totalgens = min_totalgens
	elif totalgens > max_totalgens:
		totalgens = max_totalgens
	return totalgens

def validate_irreversible_generations(irreversible_gens):
	" Checks that the value of irreversible generations is between 1000 and 500000"
	min_irreversible_gens = 1000
	max_irreversible_gens = 500000
	if irreversible_gens < min_irreversible_gens:
		irreversible_gens = min_irreversible_gens
	elif irreversible_gens > max_irreversible_gens:
		irreversible_gens = max_irreversible_gens
	return irreversible_gens

def validate_twoway_generations(twoway_gens):
	" Checks that the value of two way generations is between 1000 and 50000"
	min_twoway_gens = 1000
	max_twoway_gens = 50000
	if twoway_gens < min_twoway_gens:
		twoway_gens = min_twoway_gens
	elif twoway_gens > max_twoway_gens:
		twoway_gens = max_twoway_gens
	return twoway_gens