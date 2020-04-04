'''
Created on 5th November 2012
Code ported from Matlab code authored by Dr Matthew Hamilton of Georgetown University
General model of Natural Selection at a Dialleclic locus.
This module simulates the deterministic process of natural selection
acting on the three genotypes present at a diallelic locus.  The module
displays the frequency of the A allele over time as well and the
frequency of each genotype over time. 

User-entered parameters are:
- the initial genotype frequencies p(AA) and p(Aa), with (P(aa) determined by subtraction 
since the genotype frequencies sum to one;
- the relative fitness values of each of the three genotypes;
- the number of generations to simulate.

Output graphs show:
- Genotype frequencies over time;
- Frequency of the A allele (p) over time;
- Population mean fitness over allele frequency;
- Change in allele frequency (delta p) over allele frequency.
'''
from __future__ import division
import matplotlib.figure as fig
import numpy as np
import scipy.stats as stats
import random as rand
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

def validate_generations(generations):
     generations_min = 5
     generations_max = 5000
    
     generations = int(round(generations))
     if generations < generations_min:
        generations = generations_min
     elif generations > generations_max:
        generations = generations_max
     return generations
 
def natural_selection(generations, wAA, wAa, waa, fAA,fAa,faa):
    p = fAA+0.5*fAa
    povert = np.zeros((generations))
    wbarovert = np.zeros((generations))
    #initialize first generation
    povert[0] = p #set allele frequency for first generation
    wbarovert[0] = wAA*povert[0]**2+wAa*2*povert[0]*(1-povert[0])+waa*(1-povert[0])**2
    
    #get allele frequencies over time
    for j in xrange(1,generations):
        povert[j] = (wAA*povert[j-1]**2+wAa*povert[j-1]*(1-povert[j-1]))/wbarovert[j-1]
        wbarovert[j] = wAA*povert[j]**2+wAa*2*povert[j]*(1-povert[j])+waa*(1-povert[j])**2
    
    #get genotype frequencies over time for one initial allele frequency
    genotypeFrequencies = np.zeros((3,generations))
    for j in xrange(generations):
        genotypeFrequencies[0,j] = povert[j]**2
        genotypeFrequencies[1,j] = 2*povert[j]*(1-povert[j])
        genotypeFrequencies[2,j] = (1-povert[j])**2
    
    #get delta p over p
    points = 100
    steps = 1/points
    pt = np.arange(0,1.0,steps) #added steps because loop doesn't include end-point (want 1.0 in list)
    deltap = np.zeros((points))
    wbaroverp = np.zeros((points))
    for i in xrange(points):
        wbar = wAA*pt[i]**2+wAa*2*pt[i]*(1-pt[i])+waa*(1-pt[i])**2
        wbaroverp[i] = wbar
        deltap[i] = ((wAA*pt[i]+wAa*(1-pt[i]))/wbar)*pt[i]-pt[i]
    return (povert, genotypeFrequencies, wbaroverp,deltap,pt)
    
def plot(fAA,fAa,faa,wAA,wAa,waa,generations, clear):
    povert, genotypeFrequencies, wbaroverp,deltap,pt = natural_selection(generations, wAA, wAa, waa, fAA, fAa, faa)
    figure = fig.Figure((12, 10))
    figure.set_facecolor('#ececec')
    axes0 = figure.add_subplot(221)
    axes2 = figure.add_subplot(222)
    axes1 = figure.add_subplot(223)
    axes3 = figure.add_subplot(224)
    figure.suptitle('Natural Selection at One Locus with Two Alleles', size=16,color='#003E7E',family='serif', weight='bold')
    figure.subplots_adjust(wspace = 0.3, hspace = 0.2)
    
    if not clear:
        axes0.plot(genotypeFrequencies[0,:], 'k-.', label='AA')
        axes0.plot(genotypeFrequencies[1,:], 'b', label='Aa')
        axes0.plot(genotypeFrequencies[2,:], 'r--',label='aa')
        axes0.legend(loc='best',prop={'size':12})
        #axes0.set_xlim((1, generations+1))
        axes0.set_ylim((0.0,1.05))
        axes0.set_xlabel('Generations', fontdict={'size':14,'color':'#002147','family':'serif'})    
        axes0.set_ylabel('Genotype Frequencies', fontdict={'size':14,'color':'#002147','family':'serif'})
    
        
        axes1.plot(povert,'b')
        #axes1.set_xlim((1, generations+1))
        axes1.set_ylim((0.0,1.05))
        axes1.set_xlabel('Generations', fontdict={'size':14,'color':'#002147','family':'serif'})    
        axes1.set_ylabel('Frequency of the A allele'+r'$\ \mathit{p}$', fontdict={'size':14,'color':'#002147','family':'serif'})
    
        axes2.plot(pt,wbaroverp,'b')
        axes2.set_xlabel('Frequency of the A allele'+r'$\ \mathit{p}$', fontdict={'size':14,'color':'#002147','family':'serif'})    
        axes2.set_ylabel('Average fitness of the population', fontdict={'size':14,'color':'#002147','family':'serif'})
    
        axes3.plot(pt,deltap,'b')
         # not doing anything commented out
        #axes3.plot(pt, np.zeros((np.size(pt))),'c--')
        axes3.set_xlabel('Frequency of the A allele'+r'$\ \mathit{p}$', fontdict={'size':14,'color':'#002147','family':'serif'})    
        axes3.set_ylabel('Change in frequency of the A allele'+r'$\ (\Delta \mathit{p})$', fontdict={'size':14,'color':'#002147','family':'serif'})

     #place figure on Agg FigureCanvas
    canvas = FigureCanvas(figure)
 
    return canvas