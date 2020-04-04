'''
Created 8th November 2012
Code ported from Matlab code authored by Dr Matthew Hamilton of Georgetown University
Module to evolve a polygenic quantitative trait in a single finite population.
 Key parameters in the simulation are:
 Effective population size (Ne) that determines strength of random
 genetic drift.
 
 Mutation rate at each QTL (mu).

 Initial allele frequency at all loci (p).

 The proportion of the phenotypic variation (VP) that is caused by
 environmental variation (VE)
 
 The maximum genotypic value is the largest value the phenotype can
 take based on genotype of an individual (VE can make the value
 larger still).

 The proportion of the phenotypic distribution that does not reproduce
 sets the strength of natural selection. This is a truncation
 selection model (portion of the lower part of the phenotypic distribution
 does not reproduce) that is always directional for higher phenotypic values. 

 The number of quantitative trait loci (QTL) that contribute genetic 
 variance to trait phenotypic effects of each locus. Choices are 10,
 50 or 200 QTL. As the number of QTL increase each QTL has a smaller
 amount of influence on the phenotype. 
 
 There are three options for QTL architecture. All QTL can have equal
 effects on phenotype, the effects of QTL can be gamma distributed,
 or there can be a number of QTL with major effects (the rest have
 equal effects that explain the remaining genotypic variance) where
 the effect of EACH major QTL is specified by user input.
 
 The number of generations for the simulation.
 
 The number of randomly sampled QTL to display in the first plot of
 allele frequency by generations.

 All loci are diallelc. 
%
'''
from __future__ import division #using python 3.0 behaviour of division which will always return a float instead of casting to an int
import matplotlib.pyplot as plt
import matplotlib.figure as fig
import numpy as np
import random as rand
from math import sqrt, floor, ceil
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

def plot(populationSize, numberOfQuantitativeTraitLoci, mutationRate, generationsToSimulate, 
					   initialAlleleFrequency, environmentalVariation, maxGenotypicValue, selectionThreshold, 
					   qtlArchitecture,geneticVariation,numberOfMajorQtl, numberOfLociInAlleleFrequency, clear):
	#print 'populationSize = %d\n, numberOfQuantitativeTraitLoci=%d\n, mutationRate=%f\n, generationsToSimulate=%d\n,initialAlleleFrequency=%f\n, environmentalVariation=%f\n, maxGenotypicValue=%d\n, selectionThreshold=%f\n, qtlArchitecture=%d\n,geneticVariation=%f\n,numberOfMajorQtl=%d\n, numberOfLociInAlleleFrequency=%d\n, clear=%s' %(populationSize, numberOfQuantitativeTraitLoci, mutationRate, generationsToSimulate, initialAlleleFrequency, environmentalVariation, maxGenotypicValue, selectionThreshold, qtlArchitecture,geneticVariation,numberOfMajorQtl, numberOfLociInAlleleFrequency, clear)
	figure = fig.Figure((13,15))
	figure.set_facecolor('#ececec')
	
	ax0 = figure.add_subplot(411)
	ax1 = figure.add_subplot(412)
	ax2 = figure.add_subplot(413)
	ax3 = figure.add_subplot(414)
	figure.suptitle('Quantitative Trait Evolution In A Finite Population', size=20,color='#003E7E',family='serif', weight='bold')
	figure.subplots_adjust(hspace = 0.2)
	#print ' clear %s' %clear
	if not clear:
		frequencyA, phenotypicMean, genotypicMean, phenotypicVariance, genotypicVariance, heritability = quantitative_trait(populationSize, numberOfQuantitativeTraitLoci, mutationRate, generationsToSimulate, 
		initialAlleleFrequency, environmentalVariation, maxGenotypicValue, selectionThreshold, qtlArchitecture,geneticVariation,numberOfMajorQtl)
#		 print 'freqA %s' %frequencyA
#		 print 'phenotypmean %s' %phenotypicMean
#		 print 'genotypmean %s' %genotypicMean
#		 print 'phenotypvari %s' %phenotypicVariance
#		 print 'genotypvari %s' %genotypicVariance
#		 print 'hertit %s' %heritability
		ax0.plot(range(generationsToSimulate),frequencyA[:,0],'k')
		ax0.set_ylabel('Frequency of the plus\nallele at (some) QTL', fontdict={'size':16,'color':'#002147','family':'serif'}) 
		ax0.set_xlim(0,generationsToSimulate-1)
		ax0.set_ylim(-0.02,1.04)
		colours = ['k','b', 'g', 'r', 'c', 'm', 'y']
		lengthOfColours = len(colours)
		for i in xrange(1,numberOfLociInAlleleFrequency):
			ax0.plot(range(generationsToSimulate),frequencyA[:,i],colours[i%lengthOfColours])
		
			
		ax1.plot(range(generationsToSimulate),phenotypicMean[:,0],'k', label='Mean phenotypic value')
		ax1.plot(range(generationsToSimulate),genotypicMean[:,0],'b',label='Mean genotypic value')
		ax1.set_ylabel('Mean values', fontdict={'size':16,'color':'#002147','family':'serif'})	
		ax1.set_xlim(0,generationsToSimulate-1)
		ax1.set_ylim(0.0,maxGenotypicValue+0.1)
		ax1.legend(loc='best',prop={'size':14})
		
		ax2.plot(range(generationsToSimulate),phenotypicVariance[:,0],'k', label='Phenotypic variance'+r'$\ \mathit{VP}$')
		ax2.plot(range(generationsToSimulate),genotypicVariance[:,0],'b', label='Genotypic variance'+r'$\ \mathit{VG}$')
		ax2.set_ylabel('Variance in value\namong individuals', fontdict={'size':16,'color':'#002147','family':'serif'}) 
		ax2.legend(loc='best',prop={'size':14})
		
		ax3.plot(range(generationsToSimulate),heritability[:,0],'k')
		ax3.set_ylabel('Heritability '+r'$\frac{VG}{VP}$', fontdict={'size':16,'color':'#002147','family':'serif'}) 
		ax3.set_xlim(0,generationsToSimulate)
		ax3.set_ylim(-0.02,1.04)
		#figure.tight_layout() 
	  
	 #place figure on Agg FigureCanvas
	canvas = FigureCanvas(figure)
 
	return canvas
def quantitative_trait(populationSize, numberOfQuantitativeTraitLoci, mutationRate, generationsToSimulate, 
					   initialAlleleFrequency, environmentalVariation, maxGenotypicValue, selectionThreshold, 
					   qtlArchitecture,geneticVariation,numberOfMajorQtl):
	#rand stuff
	
	#set phenotypic effects of QTL
	alphas = qtl_effects(qtlArchitecture, numberOfMajorQtl, geneticVariation,numberOfQuantitativeTraitLoci, maxGenotypicValue)
	#genotypes = np.zeros((populationSize,numberOfQuantitativeTraitLoci))
	genotypes = initialize_genotypes(populationSize, numberOfQuantitativeTraitLoci, initialAlleleFrequency)
	frequencyA = np.zeros((generationsToSimulate, numberOfQuantitativeTraitLoci))
	frequencyA[0,:] = get_Allele_Frequencies(genotypes, populationSize,numberOfQuantitativeTraitLoci)
	
	individualGenotypicValue = np.zeros((populationSize, generationsToSimulate))
	individualPhenotypicValue = np.zeros((populationSize,generationsToSimulate))
	phenotypicMean = np.zeros((generationsToSimulate,1)) #each row is an array, row vector. could do (x) but any references must have same shape
	genotypicMean = np.zeros((generationsToSimulate,1))
	phenotypicVariance = np.zeros((generationsToSimulate,1))
	genotypicVariance = np.zeros((generationsToSimulate,1))

	#for colulm vector can access using index [i] or index and 0 axis [i,0]
	individualGenotypicValue[:,0], individualPhenotypicValue[:,0], phenotypicMean[0,0],genotypicMean[0,0],phenotypicVariance[0,0],genotypicVariance[0,0] = get_Means_and_Variances(alphas, genotypes, numberOfQuantitativeTraitLoci, populationSize, environmentalVariation)
	
	for i in xrange(1,generationsToSimulate):
		sortedIndividualPhenotypicValue, sortedGenotypes = sort_phenotypic_value(individualPhenotypicValue, genotypes, populationSize, i-1, numberOfQuantitativeTraitLoci)
	
		genotypes = selection_and_reproduction(sortedGenotypes, populationSize, numberOfQuantitativeTraitLoci, selectionThreshold, i-1)
	   # print' loop geno %s' %genotypes
		
		genotypes = mutation(mutationRate, genotypes, populationSize, numberOfQuantitativeTraitLoci)
	   # print' loop 2geno %s' %genotypes
		individualGenotypicValue[:,i], individualPhenotypicValue[:,i], phenotypicMean[i,0],genotypicMean[i,0],phenotypicVariance[i,0],genotypicVariance[i,0] = get_Means_and_Variances(alphas, genotypes, numberOfQuantitativeTraitLoci, populationSize, environmentalVariation)
		
		frequencyA[i,:] = get_Allele_Frequencies(genotypes, populationSize,numberOfQuantitativeTraitLoci)
		
	
	heritability = np.zeros((generationsToSimulate,1)) # (gens,1)
	for i in xrange(generationsToSimulate):
	   # print 'genotypicVariace[%i] = %f' %(i,genotypicVariance[i,0])
		#print 'phenotypicVariace[%i] = %f' %(i,phenotypicVariance[i,0])
		heritability[i,0] = genotypicVariance[i,0]/phenotypicVariance[i,0]
	
	return (frequencyA, phenotypicMean, genotypicMean, phenotypicVariance, genotypicVariance, heritability)

def mutation(mutationRate, genotypes, populationSize, numberOfQuantitativeTraitLoci):
	newGenotypes = genotypes
	expMuts = mutationRate*(2*populationSize)*numberOfQuantitativeTraitLoci
	numberOfMutations = int(floor(expMuts))
	fractionalNumMuts = expMuts % 1
	chance = rand.random()
	
	if fractionalNumMuts >= chance:
		numberOfMutations += 1
	
	#print ' shape of genotypes %s' %(np.shape(genotypes),)
	#print ' shape of new genotypes %s' %(np.shape(newGenotypes),)
	#choose random loci within random individuals to mutate
	for i in xrange(numberOfMutations):
		#print 'i=%d' %i
		individual = int(floor(populationSize*rand.random()))
		locus = int(floor(numberOfQuantitativeTraitLoci*rand.random()))
	   # print 'individual = %d' %individual
		#print 'locus = %d' %locus
		if genotypes[individual,locus] == 1: #genotype experiencing mutation is a 1
			newGenotypes[individual, locus] = 2 #change to 2 or Aa genotype since one gamete mutated
		elif genotypes[individual, locus] == 2: #change to 1 (AA) or 3 (aa) with 50% propbability since one gamete mutated
			chance = rand.random()
			if chance >= 0 and chance <= 0.5:
				newGenotypes[individual,locus] = 1
			else:
				newGenotypes[individual,locus] = 3 
		elif genotypes[individual,locus] == 3: # change to 2 or Aa genotype since one gamete mutated
			newGenotypes[individual,locus] = 2
   
	return newGenotypes

def selection_and_reproduction(sortedGenotypes, populationSize, numberOfQuantitativeTraitLoci, selectionThreshold,generationsToSimulate):
	'''
	make progeny based on natural selection against individuals at bottom
	of sorted list phenotypic values.  Get genotypes of new generation.
	Genotype codes:
	 1 = AA 
	 2 = Aa
	 3 = aa
	'''
	numberOfParents = populationSize - round(populationSize*selectionThreshold)
	newGenotypes = np.zeros((populationSize,numberOfQuantitativeTraitLoci))
	
	for i in xrange(populationSize):
		parent1 = int(ceil(numberOfParents*rand.random()))
		parent2 = int(ceil(numberOfParents*rand.random()))
		
		while parent1 == parent2:
			parent2 = int(floor(populationSize*rand.random())) #floor rather than ceil, need to stay below population size
		
		for j in xrange(numberOfQuantitativeTraitLoci):
#			 print ' sortedgp1[%d,%d]' %(parent1,j)
#			 print ' sortedgp2[%d,%d]' %(parent2,j)
			if sortedGenotypes[parent1, j] == 1 and sortedGenotypes[parent2,j] ==1:
				newGenotypes[i,j] = 1
			elif sortedGenotypes[parent1, j] == 2 and sortedGenotypes[parent2,j] ==2:
				#expected progeny genotypes are 1/4 1, 1/2 2 and 1/4 3
				chance = rand.random()
				if chance >= 0 and chance <= 0.25:
					newGenotypes[i,j] = 1
				elif chance > 0.25 and chance <= 0.75:
					newGenotypes[i,j] = 2
				elif chance > 0.75:
					newGenotypes[i,j] = 3
			elif sortedGenotypes[parent1, j] == 3 and sortedGenotypes[parent2,j] ==3:
				#progeny genotype is 3
				newGenotypes[i,j] = 3
			elif (sortedGenotypes[parent1, j] == 1 and sortedGenotypes[parent2,j] ==2) or (sortedGenotypes[parent1, j] == 2 and sortedGenotypes[parent2,j] ==1):
				#expected progeny genotypes are 1/2 1, 1/2 2
				chance = rand.random()
				if chance >= 0 and chance <= 0.5:
					newGenotypes[i,j] =1
				else:
					newGenotypes[i,j] = 2
			elif (sortedGenotypes[parent1, j] == 1 and sortedGenotypes[parent2,j] ==3) or (sortedGenotypes[parent1, j] == 3 and sortedGenotypes[parent2,j] ==1):
				newGenotypes[i,j] = 2
			elif (sortedGenotypes[parent1, j] == 2 and sortedGenotypes[parent2,j] ==3) or (sortedGenotypes[parent1, j] == 3 and sortedGenotypes[parent2,j] ==2):
				#expected progeny genotypes are 1/2 1, 1/2 2
				chance = rand.random()
				if chance >= 0 and chance <= 0.5:
					newGenotypes[i,j] =2
				else:
					newGenotypes[i,j] = 3
	   
	return newGenotypes
			
			
def sort_phenotypic_value(individualPhenotypicValue, genotypes, populationSize, generationToSimulate, numberOfQuantitativeTraitLoci):	
	'''
	 Make sorted array of individuals based on their phenotypic values. Phenotypic values are sorted into descending array.
	 Then make ordered list of individual genotypes based on sort order of phenotypes
	 genotypes
	 individualphenotypicvalue[:,generationToSimulate]
	
	[sortedindividualphenotypicvalue,sindices] = sort(individualphenotypicvalue[:,gen],'descend');	 sort individuals by their phenotypic values
	
	sortedgenotypes = genotypes(sindices,:);   use sorted index of phenotypic values to sort genotypes the same way

	'''
#	 print 'genotypes- sort_phenotypiv_value %s' %genotypes
#	 print 'individ[0,%d]' %generationToSimulate
#	 #find largest phenotypic value
	refValue = individualPhenotypicValue[0,generationToSimulate] #first phenotypic value
#	 print 'refvalue %s' %refValue
	refIndex = np.zeros((populationSize,1), dtype=int)
	refIndex[0,0] =0 #store index of refvalue
	
	for i in xrange(populationSize):
		if individualPhenotypicValue[i,generationToSimulate] > refValue:
			refValue = individualPhenotypicValue[i,generationToSimulate] #assign new greater phenotypic value
			refIndex[0,0] = i # store row index	 of new refvalue
	
	sortedIndividualPhenotypicValue = np.zeros((populationSize,1))
	sortedIndividualPhenotypicValue[0,0] = refValue # assign the largest value to the first entry in sorted list
	
	#find smallest	phenotypic value
	refValue = individualPhenotypicValue[0, generationToSimulate]
	refIndex[populationSize-1,0] =1
	
	for i in xrange(populationSize):
		if individualPhenotypicValue[i,generationToSimulate] < refValue:
			refValue = individualPhenotypicValue[i,generationToSimulate] #assign new smaller phenotypic value
			refIndex[populationSize-1,0] = i # store row index	of new refvalue
			
	sortedIndividualPhenotypicValue[populationSize-1,0]=refValue # assign the smallest value to the last entry in sorted list
	
	#sort phenotypic values between the largest and smallest
	for i in xrange(populationSize-2):
		startDiff = sortedIndividualPhenotypicValue[i,0] - sortedIndividualPhenotypicValue[populationSize-1,0]
		refValue = individualPhenotypicValue[i,generationToSimulate] # assign new phenotypic value
		refIndex[i+1] = i # store row index of new refvalue
		
		for j in xrange(populationSize):
			difference = sortedIndividualPhenotypicValue[i,0] - individualPhenotypicValue[j,generationToSimulate]
			if difference < startDiff and difference > 0:
				startDiff = difference #update reference difference to new small value
				refValue = individualPhenotypicValue[j,generationToSimulate]
				refIndex[i+1,0] = j
		sortedIndividualPhenotypicValue[i+1] = refValue
	
	sortedGenotypes = genotypes
	for i in xrange(populationSize):
		row = refIndex[i,0]
		
		for j in xrange(numberOfQuantitativeTraitLoci):
			sortedGenotypes[i,j] = genotypes[row,j]
			
	return (sortedIndividualPhenotypicValue, sortedGenotypes)

def get_Means_and_Variances(alphas, genotypes, numberOfQuantitativeTraitLoci, populationSize, environmentalVariation):
	'''
	 Calculate geneotypic means and variances from QTL genotypes.
	 Locus effects are vector alphas.
	 The mean phenotype for the "plus allele" homozygous genotype at locus i is equal to alpha(1,i) 
	 Heterozygotes are at the midpoint (1/2 alpha) and the "minus allele" homozygotes have a value of zero.	 
	 This is equivalent to assuming strictly additive gene action (no dominance, no epistasis). 
	'''
	 #Genotypic values are a function of locus effects alone.
	 #Phneotpic values are locus effects plus the random deviation of
	 #environmental variation.
	genotypicValue = np.zeros((populationSize, numberOfQuantitativeTraitLoci))
	phenotypicValue =  np.zeros((populationSize, numberOfQuantitativeTraitLoci))
	for i in xrange(numberOfQuantitativeTraitLoci):
		for j in xrange(populationSize):
			environmentalDeviation = np.random.randn()* sqrt(environmentalVariation)*alphas[i]
			
			if genotypes[j,i] ==1: # genotype AA
				genotypicValue[j,i] = alphas[i]
				phenotypicValue[j,i] = alphas[i] + environmentalDeviation
			elif genotypes[j,i] == 2:
				genotypicValue[j,i] = 0.5*alphas[i]
				phenotypicValue[j,i] = 0.5 * alphas[i] + environmentalDeviation
			else:
				genotypicValue[j,i] = 0.0
				phenotypicValue[j,i] = environmentalDeviation
	
	individualGenotypicValues = np.zeros((populationSize,1))
	individualPhenotypicValues	= np.zeros((populationSize,1))
	
	totalSumOfGenotypicValue = 0.0
	totalSumOfPhenotypicValue = 0.0
	
	for i in xrange(populationSize):
		lociSumOfGenotypicValue = 0.0
		lociSumOfPhenotypicValue = 0.0
		
		for j in xrange(numberOfQuantitativeTraitLoci):
			lociSumOfGenotypicValue += genotypicValue[i,j] 
			lociSumOfPhenotypicValue += phenotypicValue[i,j]
		
		individualGenotypicValues[i,0] = lociSumOfGenotypicValue
		individualPhenotypicValues[i,0] = lociSumOfPhenotypicValue
		totalSumOfGenotypicValue +=individualGenotypicValues[i,0]
		totalSumOfPhenotypicValue += individualPhenotypicValues[i,0] 
	
   # print 'total summ G %f' %totalSumOfGenotypicValue
   # print 'GI matrix %s' %individualGenotypicValues
	#print 'total summ P %f' %totalSumOfPhenotypicValue
	genotypicMean = totalSumOfGenotypicValue / populationSize
	phenotypicMean = totalSumOfPhenotypicValue / populationSize
	#print 'gmean %f' %genotypicMean
	phenotypicVariance = np.var(individualPhenotypicValues)
	genotypicVariance = np.var(individualGenotypicValues)	
	#print 'pvariane %f' %phenotypicVariance
	#print 'gvariance %f' %genotypicVariance	   
	return (individualGenotypicValues[:,0], individualPhenotypicValues[:,0],phenotypicMean,genotypicMean, phenotypicVariance, genotypicVariance)
   # other option would be to define
   #individualGenotypicValues = np.zeros((populationSize))
	#individualPhenotypicValues	 = np.zeros((populationSize))
	#and return
   #return (individualGenotypicValues,individualPhenotypicValues,phenotypicMean,genotypicMean,phenotypicVariance,genotypicVariance)
def initialize_genotypes(populationSize, numberOfQualitativeTraitLoci, initialAlleleFrequency):
	'''
	Specify the initial genotype frequencies. 
	Initial genotypes are in H-W expected frequencies based on an allele frequency p
	Genotype codes:
	1 = AA 
	2 = Aa
	3 = aa 
	'''	  
	order = np.zeros((2,3), dtype=int)
	expectedAAFrequency = int(round(initialAlleleFrequency**2*populationSize)) 
	expectedAaFrequency = int(round(2*(initialAlleleFrequency*(1-initialAlleleFrequency))*populationSize))
	expectedaaFrequency = int(round((1-initialAlleleFrequency)**2*populationSize))
#	 print 'expectedAAFrequency %d' %expectedAAFrequency 
#	 print 'expectedAaFrequency %d' %expectedAaFrequency 
#	 print 'expectedaaFrequency %d' %expectedaaFrequency 
	diff = populationSize - (expectedAAFrequency+expectedAaFrequency+expectedaaFrequency)
	if diff != 0:
		expectedaaFrequency += diff
	genotypes = np.zeros((populationSize, numberOfQualitativeTraitLoci)) 
	for i in xrange(numberOfQualitativeTraitLoci):
		# determine random order for assiging genotypes to each locus
		chance = rand.random()
		if chance >= 0 and chance < 0.3333:
			order[0,0] = expectedAAFrequency
			order[1,0] = 1
			order[0,1] = expectedAaFrequency
			order[1,1] = 2
			order[0,2] = expectedaaFrequency
			order[1,2] = 3
		elif chance >= 0.3333 and chance < 0.6667:
			order[0,1] = expectedAAFrequency
			order[1,1] = 1
			order[0,2] = expectedAaFrequency
			order[1,2] = 2
			order[0,0] = expectedaaFrequency
			order[1,0] = 3
		else:
			order[0,2] = expectedAAFrequency
			order[1,2] = 1
			order[0,0] = expectedAaFrequency
			order[1,0] = 2
			order[0,1] = expectedaaFrequency
			order[1,1] = 3
	   # print 'orders %s' %order
		counter = 0
		
		#print ' genotypes.shape %s' %(np.shape(genotypes),)
		for j in xrange(order[0,0]):
		   # print 'rows: genotypes[%d,%d] = %d' %(j,i,order[1,0])
			genotypes[j,i] = order[1,0]
			counter += 1
			
		for j in xrange( counter, counter+order[0,1]):
			#print ' counter: genotypes[%d,%d] = %d' %(j,i,order[1,1])
			genotypes[j,i] = order[1,1]
			counter += 1
		for j in xrange(counter, populationSize):
			#print 'counter2: genotypes[%d,%d] = %d' %(j,i,order[1,2])
			genotypes[j,i] = order[1,2]
	return genotypes

def get_Allele_Frequencies(genotypes, populationSize, numberOfQuanlitativeTraitLoci):
	'''
	Calculate allele frequencies from genotype matrix
	'''
	frequency_A_Allele = np.zeros((numberOfQuanlitativeTraitLoci))
	for i in xrange(numberOfQuanlitativeTraitLoci):
		countA = 0
		for j in xrange(populationSize):
			if genotypes[j,i] == 1:
				countA += 2
			elif genotypes[j,i] == 2:
				countA += 1
		frequency_A_Allele[i] = countA / (2*populationSize)
	return frequency_A_Allele	 
		
def qtl_effects(qtlArchitecture, numberOfMajorQtl, geneticVariation,numberOfQuantitativeTraitLoci, maxGenotypicValue):
	'''
	 Specify the distribution of QTL effects in the phenotype. 
	 There are three fixed QTL effect distributions:
	   effect_set = 1: equal effects at all loci

	   effect_set = 2: nMajorLociloci with large effects (the proportion of VG is specified by major_locus_value) and
	   the rest of the loci with equal effects

	   effect_set = 3: gamma distributed effect sizes
	'''			
	alphas = np.zeros((numberOfQuantitativeTraitLoci))
	if qtlArchitecture == 1:
		for i in xrange(numberOfQuantitativeTraitLoci):
			alphas[i] = maxGenotypicValue/numberOfQuantitativeTraitLoci
	elif qtlArchitecture == 2:
		propOfEquals = 1 - numberOfMajorQtl * geneticVariation
		effect = (maxGenotypicValue * propOfEquals)/(numberOfQuantitativeTraitLoci-numberOfMajorQtl)
	   # print 'alphas : %s' %alphas
		for i in xrange(numberOfMajorQtl):
#			 print 'alphas : %s' %alphas
#			 print 'i=%d' %i
#			 print 'alphas[%s] = %s' %(i,alphas[i])
#			 print 'geneticvariation * maxgenotypevlaue -> %s * %s = %s' %(geneticVariation,maxGenotypicValue,geneticVariation * maxGenotypicValue)
			alphas[i] = geneticVariation * maxGenotypicValue
			
		for i in xrange(numberOfMajorQtl, numberOfQuantitativeTraitLoci):
			alphas[i] = effect
	elif qtlArchitecture == 3:
		#shape parameters for gamma distribution
		gammaA = 1
		gammaB = 1
		alphas = np.random.gamma(gammaA,gammaB,(numberOfQuantitativeTraitLoci))
		sum = 0.0
		for i in xrange( numberOfQuantitativeTraitLoci):
			sum += alphas[i]
		scale = maxGenotypicValue/sum
		
		for i in xrange(numberOfQuantitativeTraitLoci):
			alphas[i] = alphas[i] * scale
	
	return alphas

def validate_population_size(populationSize):
	" Checks that the value of populationSize is between 1 and 10000"
	print 'incoming population size %d' %populationSize
	min_populationSize = 1
	max_populationSize = 10000
	populationSize = int(round(populationSize))
	if populationSize < min_populationSize:
		populationSize = min_populationSize
	elif populationSize > max_populationSize:
		populationSize = max_populationSize
	print 'outgoing population size %d' %populationSize
	return populationSize

def	 validate_generations(generations):
	print 'incoming generations %d' %generations
	generations_min = 10
	generations_max = 10000
	generations = int(round(generations))
	if generations < generations_min:
		generations = generations_min
	elif generations > generations_max:
		generations = generations_max
		
	print 'outgoing generations %d' %generations
	return generations

def validate_number_of_major_qtl(qtlNumber):
	min_qtlNumber = 1
	max_qtlNumber = 6
	
	qtlNumber = int(round(qtlNumber))
	if qtlNumber < min_qtlNumber:
		qtlNumber = min_qtlNumber
	elif qtlNumber > max_qtlNumber:
		qtlNumber = max_qtlNumber
	return	qtlNumber

def validate_genotypic_variance(genotypicVariance, qtlNumber):
	min_genotypicVariance = 0
	max_genotypicVariance = 1/qtlNumber
	
	if genotypicVariance < min_genotypicVariance:
		genotypicVariance = min_genotypicVariance
	elif genotypicVariance > max_genotypicVariance:
		genotypicVariance = max_genotypicVariance
	
	return genotypicVariance
	