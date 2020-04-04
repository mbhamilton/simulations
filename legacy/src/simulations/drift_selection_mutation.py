'''
Created on  22 October, 2012
Code ported from Matlab code authored by Dr Matthew Hamilton of Georgetown University
The code below draws the H-W parabola on a definetti diagram and plots points that represent
a set of genotype frequencies for one locus with two alleles
@author: marieidowu
'''
from __future__ import division #using python 3.0 behaviour of division which will always return a float instead of casting to an int
import matplotlib.figure as fig
import numpy as np
import random as rand
from math import sqrt
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import matplotlib.gridspec as gridspec 

def genetic_drift(pmatrix, ne, j, nreps):
    #print pmatrix
    for i in xrange(nreps):
       # print ('pmatrix[%d,%d]=%f') %(i,j,pmatrix[i,j])
        delta = (2*rand.random()-1)* sqrt((3*pmatrix[i,j]*(1-pmatrix[i,j]))/(2*ne))
        #print'delta %f' %delta
        if pmatrix[i,j] + delta >1.0:
            pmatrix[i,j]=1.0
        elif pmatrix[i,j] + delta < 0.0:
            pmatrix[i,j]=0.0
        else:
            pmatrix[i,j] = pmatrix[i,j]+delta
    return pmatrix

def natural_selection(pmatrix, gen, wAA,wAa,waa,nreppops):
#    print 'pmatrix %s' %pmatrix
#    print 'gen %d' %gen
#    print 'wAA %s' %wAA
#    print 'wAa %s' %wAa
#    print 'waa %s' %waa
#    print 'nreppops %s' %nreppops
    for j in xrange(nreppops):
        wbar = wAA*pmatrix[j,gen]**2+wAa*2*pmatrix[j,gen]*(1-pmatrix[j,gen])+waa*(1-pmatrix[j,gen])**2
        pmatrix[j,gen] = (wAA*(pmatrix[j,gen]**2)+wAa*pmatrix[j,gen]*(1-pmatrix[j,gen]))/wbar
        
   # print 'natrual selection method'
    #print 'natural selection method return \npmatrix = %s'  %pmatrix
    return  pmatrix


def get_number_of_mutations(eff_mut_rate):
    randomnumber = rand.random()
    if eff_mut_rate >= 1.0:
        number_of_mutations = round(eff_mut_rate)
        if (eff_mut_rate % 1 >= rand.random()):
            number_of_mutations += 1
    elif eff_mut_rate >= randomnumber:
        number_of_mutations = 1
    else:
        number_of_mutations = 0
   # print 'number of mutations %d' %number_of_mutations
    return number_of_mutations

def mutation(pmatrix,ne,gen,mu,nu,nreppops):
    #print np.shape(pmatrix)
    freq_one_allele = 1/(2*ne)
    eff_forward_mututation_rate = 2*ne*mu
    for i in xrange(nreppops):
       # print 'i = %d' %i
        #print 'pmatrix[%d,%d]' %(i,gen)
        number_of_forward_mutations = get_number_of_mutations(eff_forward_mututation_rate)
            
        if number_of_forward_mutations > 0.0:
            delta_p = number_of_forward_mutations*freq_one_allele
            if pmatrix[i,gen]-delta_p < 0.0:
                pmatrix[i,gen] = 0.0
            else:
                pmatrix[i, gen] = pmatrix[i,gen] - delta_p
                
        eff_reverse_mututaion_rate = 2*ne*nu
        number_of_reverse_mutations = get_number_of_mutations(eff_reverse_mututaion_rate)
        
        if number_of_reverse_mutations > 0.0:
            delta_p = number_of_reverse_mutations * freq_one_allele
            if pmatrix[i,gen]+delta_p > 1.0:
                pmatrix[i,gen] = 1.0
            else:
                pmatrix[i,gen] = pmatrix[i,gen]+delta_p
    return pmatrix
                
def validate_number_of_replicate_populations(number_of_populations):
    min_number_of_populations =20
    max_number_of_populations =10000
    number_of_populations = int(round(number_of_populations))
    
    if number_of_populations < min_number_of_populations:
        number_of_populations = min_number_of_populations
    elif number_of_populations > max_number_of_populations:
        number_of_populations= max_number_of_populations
        
    return number_of_populations

def validate_value(value):
    min_value = 5
    max_value = 10000
    
    value = int(round(value))
    
    if value < min_value:
        value = min_value
    elif value > max_value:
        value = max_value
    
    return value
    
def plot(number_of_replicate_populations, population_size, initial_allele_frequency, number_of_generations, 
         number_of_replicate_populations_in_allele,natural_selection_selected, wAA,wAa,waa,mutation_selected,mu,nu,clear ):
    pmatrix = np.zeros((number_of_replicate_populations,number_of_generations))
    pmatrix[:,0] = initial_allele_frequency
    for j in xrange(1, number_of_generations):
        pmatrix[:,j] = pmatrix[:,j-1]
      #  print 'pmatrix in plot metnod  %s' %pmatrix
        if mutation_selected:
            pmatrix = mutation(pmatrix, population_size,j,mu,nu,number_of_replicate_populations)
        
        pmatrix = genetic_drift(pmatrix, population_size, j, number_of_replicate_populations)
        
        if natural_selection_selected:
#            print 'selection = %s' %natural_selection_selected
            pmatrix = natural_selection(pmatrix, j, wAA, wAa, waa, number_of_replicate_populations)
    
    figure = fig.Figure((15,12))
    figure.set_facecolor('#ececec')
    figure.suptitle('Genetic drift (+ selection + mutation) in many independent replicate populations', size=18,color='#003E7E',family='serif', weight='bold')

    gs = gridspec.GridSpec(2,5)#, width_ratios=[2,1],height_ratios=[2,1])
    gs.update(wspace=0.5) # This is similar to subplots_adjust, but it only affects the subplots that are created from the given Gridspec
    axes = figure.add_subplot(gs[0:,:-2])#2]) #3x3 subplot
    axes2 = figure.add_subplot(gs[1:,3:])
    axes2.set_title('Distribution of allele frequencies\n among all replicate populations \nat last generation\n', fontdict={'size':16,'color':'#002147','family':'serif','weight':'bold'})

    if not clear:
        replacements = rand.sample(xrange(1,number_of_replicate_populations_in_allele+1),number_of_replicate_populations_in_allele)
        axes.plot(pmatrix[replacements[0],:],'b-')
        
        for i in xrange(1,number_of_replicate_populations_in_allele):
            axes.plot(pmatrix[replacements[i],:],'b-')
       
        axes.set_ylabel('Frequency of A allele (p)', fontdict={'size':14,'color':'#002147','family':'serif'})
        axes.set_xlabel('Generations', fontdict={'size':14,'color':'#002147','family':'serif'})
        #axes.set_xlim(0,number_of_generations)
        axes.set_ylim(-0.01,1.01)
        
        axes2.set_xlabel('Frequency of A allele (p)', fontdict={'size':14,'color':'#002147','family':'serif'})
        axes2.set_ylabel('Number of replicate populations', fontdict={'size':14,'color':'#002147','family':'serif'})
      
        #print 'pmatrix = %s' %(pmatrix,)
        #print 'pmatrix[:,%d] = %s' %(number_of_generations-1,pmatrix[:,number_of_generations-1])
        axes2.hist(pmatrix[:,number_of_generations-1])#, bins=np.arange(0,1, 0.1))
        #axes2.set_xlim(0.0,1)
    canvas = FigureCanvas(figure)

    return canvas
   
