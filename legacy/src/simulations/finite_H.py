'''
Created on October 01, 2012
Code ported from Matlab code authored by Dr Matthew Hamilton of Georgetown University
Decline in heterosygosity  over time in a finite population
 
 This simulation demonstrates how genetic drift in finite populations decreases heterozygosity over time. 
 In this sense, genetic drift and consanguineous mating have the same impact on heterozygosity. 
 
 User-entered parameters are the 'Effective population size (Ne)', the 'Initial allele frequency p(A)',
 the 'Number of generations' to simulate, and the 'Number of replicate populations to plot.' 
 
 Press Run to generate a graph of heterozygosity (H) over time in generations. 
 
 In a finite population, expected heterozygosity decreases according to
 Ht = H0(1-1/(2Ne))^t
 where Ht is heterozygosity at time t, H0 is the initial heterozygosity at time zero, t is the number
 of generations, and Ne is the effective population size. The generation-by-generation values computed
 from this equation are shown as a solid blue line on the graph. 
 
 There are additional lines on the graph (the number of lines is equal to the value in the text entry box 
 ?Number of replicate populations to plot?). These lines each show the heterozygosity over time in a single 
 finite population that has an effective population size equal to the value in the text entry field. Notice 
 that the heterozygosity frequencies in the individual populations varies over time and is often greater or 
 less than the expected heterozygosity. This reinforces that the equation for expected heterozygosity applies 
 to the average of a very large number of replicate populations (it is an ensemble distribution). Individual 
 populations do experience a reduction in heterozygosity over time caused by genetic drift, but the decline 
 in heterozygosity is a random walk. 
@author: marieidowu
'''

from __future__ import division #using python 3.0 behaviour of division which will always return a float instead of casting to an int
import matplotlib.figure as fig
import numpy as np
import random as rand
from math import sqrt, floor
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import logging as lg


def drift_h(  population_size, number_of_generations, number_of_generations_between_mutations, number_of_replicate_populations, clear):
    figure = fig.Figure((12,10))
    figure.set_facecolor('#ececec')
    axes = figure.add_subplot(1,1,1)
    #clear_figure(figure)
    #axes = figure.gca()
    
   
    figure.suptitle('Heterozygosity at Diallelic locus in a Finte Population', size=18,color='#003E7E',family='serif', weight='bold')
    
    if not clear:
        h, exph = h_over_time(population_size, number_of_generations, number_of_generations_between_mutations, number_of_replicate_populations)
        axes.plot(exph.transpose(),'b') #exph is a 1xn matrix/array, had to transpose to get column vector for graph to plot, exph was giving a blank plot
        
        h=h.transpose()
        colours = ['c','k','m','r','g']
        symbols = ['o','x','s','d','*','v','p']
        length_of_colours = len(colours)
        length_of_symbols = len(symbols)
        #print 'h size = ',h.shape
    #    
    #    # to avoid list out of range error checking the size of lists against m and starting at the
    #    # beginning of the list once m is equal or greater than the size of the list
        for m in xrange(number_of_replicate_populations):               
    #        print 'm=%d' %m
            sym = symbols[m%length_of_symbols]
     #       print 'sym[%d]=%s' %(m%length_of_symbols,sym)
            col = colours[m%length_of_colours]
      #      print 'colors[%d]=%s' %(m%length_of_colours,col)
            linemark = ':'+col
       #     print 'linemark= %s' %linemark
            linemark = linemark+sym
        #    print 'linemark= %s' %linemark
    #        
    #        print 'ploting h[:%d], linemark %s' %(m,linemark)
    #        print 'h =', h[:,m], '\n h size =', (h[:,m]).shape
            axes.plot(h[:,m],linemark)
        axes.set_xlabel('Generation', fontdict={'size':16,'color':'#002147','family':'serif'})
        axes.set_ylabel('Heterozygote genotype frequency (pAa)', fontdict={'size':16,'color':'#002147','family':'serif'})
        axes.set_xlim(0,number_of_generations)
        axes.set_ylim(0.0,2*number_of_generations_between_mutations*(1-number_of_generations_between_mutations))
        
        #set a legend for the plot using the labels specified in the above
        #plot commands, set it in the best location with fontsize of 12
        axes.legend(loc='best',prop={'size':12})
    #place figure on Agg FigureCanvas
    canvas = FigureCanvas(figure)

    return canvas
        
def h_over_time(population_size, number_of_generations, number_of_generations_between_mutations, number_of_replicate_populations):
    '''
        Show expected decline in heterozygosity with time Ht=(1-1/(2Ne))**t*H0 (Ne=population_size) 
        compared to observe heterozygosity in population experiencing genetic drift
    '''
    q = 1-number_of_generations_between_mutations
    pmatrix = np.zeros((number_of_replicate_populations,number_of_generations))
    #initialize allele data matrix
    pmatrix[:,0] = number_of_generations_between_mutations
    
    
    #drift allele frequencies
    for i in xrange(number_of_replicate_populations):
        pmatrix = drift_allele_frequency(pmatrix, number_of_generations, population_size, i)
    #express effective population size as chnage in h over time for these same populations
    initial_h = 2*number_of_generations_between_mutations*q
    h = 2*pmatrix*(1-pmatrix) #element wise multiplication
    
    #allocate space
    exph = np.zeros((1,number_of_generations))
    exph[0,0] = initial_h
    
    for i in xrange(1,number_of_generations):
        exph[0,i] = exph[0,i-1]*(1-(1/(2*population_size)))
    
    return (h, exph)
        
def drift_allele_frequency(data_matrix, number_of_generations, population_size,row): 
    '''
        allele frequencies experience drift using Kimura's psuedo random variable
        method
    ''' 
    for i in xrange(1,number_of_generations):
        delta = (2*rand.random()-1)*sqrt((3*data_matrix[row,i-1]*(1-data_matrix[row,i-1]))/(2*population_size))
        if(data_matrix[row,i-1]+delta) > 1.0:
            data_matrix[row,i] = 1.0
        elif (data_matrix[row,i-1]+delta) < 0.0:
            data_matrix[row,i] = 0.0
        else:
            data_matrix[row,i] = data_matrix[row,i-1]+delta
            
    #I don't think its necessary to return the matrix, should be passby reference
    return data_matrix

def validate_population_size(population_size):
    min = 5
    max =10000
    
    rounded_size = int(round(population_size))
    if rounded_size < min:
        rounded_size = min
    elif rounded_size > max:
        rounded_size = max
    
    return rounded_size
    
    
def validate_number_of_generations(number_of_generations):
    min = 5
    max =5000
    
    rounded_number = int(round(number_of_generations))
    if rounded_number < min:
        rounded_number = min
    elif rounded_number > max:
        rounded_number = max
    
    return rounded_number  
    
def clear_figure(figure):
    figure.clf()    