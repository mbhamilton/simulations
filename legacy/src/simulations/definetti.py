'''
Created on May 4, 2012
Code ported from Matlab code authored by Dr Matthew Hamilton of Georgetown University
The code below draws the H-W parabola on a definetti diagram and plots points that represent
a set of genotype frequencies for one locus with two alleles
@author: marieidowu
'''
import matplotlib.figure as fig
import numpy as np

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import logging as lg


def plot(showAxis, fAA, fAa, generations, matingModel, clear):
    "Plot the definetti diagram" 
    # Define parabola of hetrozygote frequency under H-W 
    points =1000.0
    steps = 1/points
           
    x = np.arange(0,1.001,steps) #added the 0.001 as range does not include the end int, without 0.001 it would stop at .999
    y = 2*x*(1-x);


    #create a Figure instance 
    figure = fig.Figure((12,10))
    figure.set_facecolor('#ececec')
    #use the Figure to create one or more Axes or Subplot instances
    axes = figure.add_subplot(1,1,1)
    axes.set_xlim((0,1))
    axes.set_ylim((0,1))
    if not clear:
        #use the Axes instance to plot of H-W parabola
        axes.plot(x,y) 
        
        #set x and y axis limits
        
        #draw in lines for triangle of permissible frequency values
        axes.plot([0.0,0.5],[0.0,1.0],linestyle='-',color='k')
        axes.plot([0.5,1.0],[1.0,0.0],linestyle='-',color='k')
        
        #draw in lines for genotype frequencies
        if showAxis:
            plot_genotype_frequency_axes(figure)
        
        xpoint, ypoint, faa = plot_single_point_of_genotype_frequencies(figure, fAA, fAa)
        simulate_over_generations(figure,xpoint,ypoint,generations,matingModel,fAA,fAa, faa)
         #set axes labels
        axes.set_ylabel("Frequency of Aa genotype", fontdict={'size':18,'color':'#002147','family':'serif'})
        axes.set_xlabel("Frequency of A allele (p)", fontdict={'size':18,'color':'#002147','family':'serif'})
        #set a legend for the plot using the labels specified in the above
        #plot commands, set it in the best location with fontsize of 12
        axes.legend(loc='best',prop={'size':18})
    
    axes.set_title('Genotype Frequencies for one locus with two alleles on a de Finetti Diagram\n\n', fontdict={'size':18,'color':'#003E7E','family':'serif'},horizontalalignment='center', weight='bold')
   #place figure on Agg FigureCanvas
    canvas = FigureCanvas(figure)
    

    return canvas

def plot_genotype_frequency_axes(figure):
    " Add the genotype frequencylines to the plot"
    currentAxes = figure.gca()
    currentAxes.plot([0.5,0.5],[1.0,0.0],linestyle='--',color='b', label='Frequency of Aa') #labels are used for the legend
    currentAxes.plot([0.0,0.76],[0.0,0.5],linestyle='-',color='k', label='Frequency of AA')
    currentAxes.plot([0.25,1.0],[0.5,0.0],linestyle=':',color='k', label='Frequency of aa') 
   

def plot_single_point_of_genotype_frequencies(figure, fAA,fAa):
    currentAxes = figure.gca()
    faa = 1 - fAA - fAa
    #plot first genotype data point
    ypoint = fAa
    # x-axis point is defined by frequency of A allele, p
    xpoint = fAA + 0.5 * ypoint #get frequency of A, that is distance along x-axis
    currentAxes.plot(xpoint, ypoint, 'rD')
    return xpoint,ypoint,faa
    
        
def simulate_over_generations(figure, xpoint,ypoint, generations, model, fAA,fAa,faa):
        
     currentAxes = figure.gca()
     
     if generations > 1:
         #array of zeros, to use matrix need to use np.matlib.zeros((gens,2))
         #can use array here since we will not be doing matrix multiplication or multiplying by a scalar
         data = np.zeros((generations,2)) #array of gens rows and 2 columns
         data[0,0] = xpoint # 1st gen xpoint - indexing starts from 0 so this is position 1st row 1st col
         data[0,1] = ypoint # 1st gen ypoint - 1st row 2nd col
         
         if model == 0: # Random mating
            p = xpoint
            q = 1-p
            D = p**2
            H = 2*p*q
            R = q**2
            #under H-W all subsequent genotype freqs are the same so only need compute 2nd
            # generation values
            data[1:generations,:1] = D+0.5*H #row 2 to generations of column 1 (column 0 from row 1), indexing is upto but not including end column
            data[1:generations,1:] = H # column 2 ( slicing from 1 to end)
            generationPoints = currentAxes.plot(data[:,:1],data[:,1:],'-ro')
            print 'gen points ',generationPoints
         elif model == 1: #Positive Assortative
            D = np.zeros((generations,1))
            H = np.zeros((generations,1))
            R = np.zeros((generations,1))
            
            D[0,0] = fAA
            H[0,0] = fAa
            R[0,0] = faa
            
            for i in xrange(1,generations):
                D[i,0] = D[i-1,0]+0.25*H[i-1,0]
                H[i,0] = 0.5*H[i-1,0]
                R[i,0] = R[i-1,0]+0.25*H[i-1,0]
                data[i,0] = D[i,0]+0.5*H[i,0] # xpoint
                data[i,1] = H[i,0] #ypoint
            generationPoints = currentAxes.plot(data[:,:1],data[:,1:],'-ro')
            print 'gen points ',generationPoints
         elif model ==2: #Negative Assortative'
            D = np.zeros((generations,1))
            H = np.zeros((generations,1))
            R = np.zeros((generations,1))
            
            D[0,0] = fAA
            H[0,0] = fAa
            R[0,0] = faa
            
            for i in xrange(1,generations):
                #recursions from Karlin and Feldman 1968
                D[i,0] = (D[i-1,0]*H[i-1,0]*(1+R[i-1,0]))/(2*(1-D[i-1,0])*(1-H[i-1,0]))
                H[i,0] = ((D[i-1,0]*R[i-1,0])/2)*((1+H[i-1,0])/((1-D[i-1,0])*(1-R[i-1,0])))+0.5
                R[i,0] = (R[i-1,0]*H[i-1,0]*(1+D[i-1,0]))/(2*(1-R[i-1,0])*(1-H[i-1,0]))
                data[i,0] = D[i,0]+0.5*H[i,0] #xpoint
                data[i,1] = H[i,0] #ypoint
            generationPoints = currentAxes.plot(data[:,:1],data[:,1:],'-ro')
            print 'gen points ',generationPoints
         print 'leaving simulate over generations'


