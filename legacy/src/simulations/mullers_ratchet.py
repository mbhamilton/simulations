'''
Created on October, 2012
Code ported from Matlab code authored by Dr Matthew Hamilton of Georgetown University
This code simulates Muller's ratchet for a population of non-recombing haploids.
Chance of survival is fitness where fitness is (1-s)^i where i is the number of loci.
This means that natural selection takes the form of viability selection.
Family size is Poisson distributed.
@author: marieidowu
'''
from __future__ import division #using python 3.0 behaviour of division which will always return a float instead of casting to an int
import matplotlib.figure as fig
import numpy as np
import scipy as sp
import scipy.stats as stats
import random as rand
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas


def validate_population_size(n):
    n_min = 10
    n_max = 10000
    
    if n < n_min:
        n = n_min
    elif n > n_max:
        n = n_max
    return n

def  validate_generations(generations):
    generations_min = 10
    generations_max = 10000
    
    if generations < generations_min:
        generations = generations_min
    elif generations > generations_max:
        generations = generations_max
    return generations

def validate_coefficient(coefficient):
    coefficient_min = 0.0000001
    coefficient_max = 0.99
    if coefficient < coefficient_min:
        coefficient = coefficient_min
    elif coefficient > coefficient_max:
        coefficient = coefficient_max    
    return coefficient

def validate_mutation(mutation):
    mutation_min = 0.00000001
    mutation_max = 0.1
    if mutation < mutation_min:
       muutation = mutation_min
    elif mutation > mutation_max:
       mutation = mutation_max   
    return mutation

def validate_loci(loci):
    loci_min = 10
    loci_max = 1000
    if loci < loci_min:
        loci = loci_min
    elif loci > loci_max:
        loci = loci_max
    return loci

def validate_mean_family_size(mean):
    meanFamilySize_min = 0.1
    meanFamilySize_max = 50
    
    if mean < meanFamilySize_min:
        mean = meanFamilySize_min
    elif mean > meanFamilySize_max:
        mean = meanFamilySize_max
    return mean

def validate_max_family_size(max):
    maxFamilySize_min = 1
    maxFamilySize_max = 100
    if max < maxFamilySize_min:
        max = maxFamilySize_min
    elif max > maxFamilySize_max:
        max = maxFamilySize_max
    return max    
    
    
def family_size_distribution(meanFamilySize, maxFamilySize): 
    familyRange = range(maxFamilySize+1)
    probabilityFamilySize = stats.poisson.pmf(familyRange,meanFamilySize)    
    sumholder = 0.0
    cumulativeProbabilityFamilySize = []
    for i in xrange(maxFamilySize+1):
        sumholder = sumholder + probabilityFamilySize[i]
        cumulativeProbabilityFamilySize.append(sumholder)
        
    #print np.shape(cumulativeProbabilityFamilySize)
    return cumulativeProbabilityFamilySize

def mullers(generations, populationSize,loci, selectionCoefficient,mutationRate,meanFamilySize,maxFamilySize ):
    population = np.zeros((generations+1,populationSize,loci))
    mutationCount = np.zeros((populationSize,1))
    mutationOverGenerations = np.zeros((generations+1,populationSize,1))
    fitness = np.ones((populationSize,1))
    progeny = np.zeros((populationSize,loci))
    progenyList = np.zeros((populationSize,1))
    cumulativeProbabilityFamilySize = family_size_distribution(meanFamilySize, maxFamilySize)
    numberOfMutations = populationSize * mutationRate
    for generation in xrange(generations):
        
        if numberOfMutations >= 1.0:
            numberOfMutations = int(round(populationSize*mutationRate))
        elif numberOfMutations >= rand.random():
            numberOfMutations=1
        else:
            numberOfMutations=0
        
        for j in xrange(numberOfMutations):
            randomIndividual = rand.randint(0,populationSize-1) # random int between 0 and populationSize-1
            randomLocus = rand.randint(0, loci-1)
            
            if population[generation,randomIndividual, randomLocus] == 0:
                population[generation,randomIndividual, randomLocus] =1
            else:
                j-=1
            sp.sum(population[generation,:,:],1, out=mutationCount) #sum over sencond dimension rows and assign to mutationcount row vector 
            totalMutations = sp.sum(mutationCount)
            if totalMutations == populationSize*loci:
                print 'All loci have a mutation, exiting loop' 
                break
        
        mutationOverGenerations[generation+1,:,:] = mutationCount
        fitness = (1-selectionCoefficient)**mutationCount
        #print 'fitness size = %s' %(np.shape(fitness),)
        progny = 0
        parent = 0
        while progny < populationSize:
#            print 'progny=%d' %progny
#            print ' parent=%d' %parent
#            print 'fitness [%d,0]=%s' %(parent,fitness[parent,0])
            if rand.random() <= fitness[parent,0]:
                prob = rand.random()
                for k in xrange(maxFamilySize+1):
                    if prob <= cumulativeProbabilityFamilySize[k]:
                        familySize = k-1
                        break
                
                #use +1 instead of -1 to avooid infinite loop when progny is 0 and maxfamsize is 1
                # progny+familySize or progny+familySize-1 will result in infinite loop
                for k in xrange(progny,progny+familySize+1):
                    if k < populationSize:
                        progenyList[k,0] =parent
                        progny+=1
                    else:
                        break
            
            if parent < populationSize-1:
                parent+=1
            else:
                parent=0
        
        for j in xrange(populationSize):
            parent = progenyList[j,0]
            for k in xrange(loci):
                progeny[j,k] = population[generation,parent,k]
        
        population[generation+1,:,:]= progeny
    return mutationOverGenerations
        
def plot(populationSize,generations,loci, s, mu, meanfs, maxfs,clear):
     "Plot the Muller's Ratchet diagram"    
     figure = fig.Figure((12,10))
     figure.set_facecolor('#ececec')
     # create subplots and assign returned objects figure and axess
     #figure, axes = plt.subplots(2, 2)
     
     axes0 = figure.add_subplot(221)
     axes1 = figure.add_subplot(222)
     axes2 = figure.add_subplot(223)
     axes3 = figure.add_subplot(224)
     figure.suptitle('Muller\'s Ratchet Simulations', size=16,color='#003E7E',family='serif', weight='bold')
     figure.subplots_adjust(wspace = 0.2, hspace = 0.2)
        
     if not clear: 
     	mutateOverGenerations = mullers(generations, populationSize, loci, s, mu, meanfs, maxfs)  
   
         #now make some graphs of the results
     
         #determine the maximum number of mutations at final time point
         #maxbin needs to be a column not row vector(array) use reshape(-1,1) to change from row to column, 
         
        maxbin = np.max(mutateOverGenerations[generations,:,:],0)+1
        bins = range(maxbin+1)
        time = 0
         
        #plot 1
        #n,bins,patches= axes0.hist(mutateOverGenerations[time,:,:],bins)
        n, bins = np.histogram(mutateOverGenerations[time,:,:],bins)
        proportion = n/populationSize
        axes0.bar(bins[0:len(proportion)], proportion)
        axes0.set_title('Generation = '+str(time), fontdict={'size':14,'color':'#002147','family':'serif'})
        axes0.set_xlim(-0.5,maxbin+0.5)
        axes0.set_ylim(0,1)
        axes0.set_ylabel('Proportion of haploid individuals', fontdict={'size':14,'color':'#002147','family':'serif'})
     
         #plot 2
        print 'plot 2: generations*0.33=%f' %round(generations*0.33)
        time = int(round(generations*0.33))
        n, bins = np.histogram(mutateOverGenerations[time,:,:],bins)
        proportion = n/populationSize
        axes1.bar(bins[0:len(proportion)], proportion)
        axes1.set_title('Generation = '+str(time), fontdict={'size':14,'color':'#002147','family':'serif'})
        axes1.set_xlim(-0.5,maxbin+0.5)
        axes1.set_ylim(0,1)
        #plot 3
        print 'plot 3: generations*0.66=%f' %round(generations*0.66)
        time = int(round(generations*0.66))
        n, bins = np.histogram(mutateOverGenerations[time,:,:],bins)
        proportion = n/populationSize
        axes2.bar(bins[0:len(proportion)], proportion)
        axes2.set_title('Generation = '+str(time), fontdict={'size':14,'color':'#002147','family':'serif'})
        axes2.set_xlim(-0.5,maxbin+0.5)
        axes2.set_ylim(0,1)
        axes2.set_ylabel('Proportion of haploid individuals', fontdict={'size':14,'color':'#002147','family':'serif'})    
        axes2.set_xlabel('Number of mutations per genome', fontdict={'size':14,'color':'#002147','family':'serif'})
    
          #plot 4
        print 'plot 4: generations=%f' %generations
        time = generations
        n, bins = np.histogram(mutateOverGenerations[time,:,:],bins)
        proportion = n/populationSize
        axes3.bar(bins[0:len(proportion)], proportion)
        axes3.set_title('Generation = '+str(time), fontdict={'size':14,'color':'#002147','family':'serif'})
        axes3.set_xlim(-0.5,maxbin+0.5)
        axes3.set_ylim(0,1)
        axes3.set_xlabel('Number of mutations per genome', fontdict={'size':14,'color':'#002147','family':'serif'})

     #place figure on Agg FigureCanvas
     canvas = FigureCanvas(figure)
 
     return canvas                    
     
         