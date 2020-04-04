'''
Created 26 September 2013
Code translated from Matlab code authored by Dr Matthew Hamilton of Georgetown University
for the ITEL pilot project.
@author: Marie Kolawole
'''

from __future__ import division
import matplotlib.figure as fig
import numpy as np
from mutation import validate_ne, validate_generations_between
import random as rand
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import matplotlib.gridspec as gridspec
from math import sqrt

def set_generations(generations, population_size):
	'''
	generations: number of generations to view
	population_size: effective population size (Ne)
	
	Set total simulation time, expected time to fixation approaches 4Ne
	'''
	# 'generations: %d' %generations
	new_generations = 4*population_size
	if new_generations < generations:
		return generations
	else:
		return new_generations
		
def validate_population_size(ne):
	return validate_ne(ne)
	
def validate_generations(generations):
	'''Checks that the value of total generations is between 100 and 500'''
	min_generations = 100
	max_generations = 500
	if generations < min_generations:
		generations = min_generations
	elif generations > max_generations:
		generations = max_generations
	return generations

def validate_generations_between_mutations(interval):
	return validate_generations_between(interval)

def validate_number_of_mutations(number_of_loci):
	min_number_of_loci = 200
	max_number_of_loci = 2000
	if number_of_loci < min_number_of_loci:
		number_of_loci = min_number_of_loci
	elif number_of_loci > max_number_of_loci:
		number_of_loci = max_number_of_loci

	return number_of_loci
	
	
def get_selection_coefficients(number_of_loci, mutation_fitness):
	selection_coefficients = np.ones((number_of_loci,1)) # column vector of selection coefficients; 1 equals mean fitness
	fitness_upper = 1+mutation_fitness/100
	fitness_lower = 1-mutation_fitness/100
	
	for i in xrange(0,number_of_loci):
		selection_coefficients[i,0] = rand.random()*(fitness_upper-fitness_lower)+fitness_lower
		
	return selection_coefficients,fitness_upper,fitness_lower
	
def get_distribution_times(number_of_loci,generations_to_simulate, neutral_frequencies, nearly_neutral_frequencies):
	'''
	Organize data into alleles that fixed and alleles that went to loss
	so that the distributions of times to fixation and loss can be determined
	'''
	neutral_generations_to_fixation = np.zeros((number_of_loci,1))
	neutral_generations_to_loss = neutral_generations_to_fixation.copy()
	
	nearly_neutral_generations_to_fixation = neutral_generations_to_fixation.copy()
	nearly_neutral_generations_to_loss = neutral_generations_to_fixation.copy()
	for x in xrange(number_of_loci): #0-num_of_loci-1
		for y in xrange(generations_to_simulate-1,-1,-1 ): # end at -1 to include 0
			if neutral_frequencies[x,y] == 1.0:
				neutral_generations_to_fixation[x,0] = y
			elif neutral_frequencies[x,y] == 0.0:
				neutral_generations_to_loss[x,0] = y
			
			if nearly_neutral_frequencies[x,y] == 1.0:
				nearly_neutral_generations_to_fixation[x,0] = y
			elif nearly_neutral_frequencies[x,y] == 0.0:
				nearly_neutral_generations_to_loss[x,0] = y
	
	return neutral_generations_to_fixation,neutral_generations_to_loss,nearly_neutral_generations_to_fixation,nearly_neutral_generations_to_loss
	
def plot(generations,population_size, number_of_loci, mutation_fitness, interval, clear):
	'''
	'''
	figure = fig.Figure((16,13),facecolor='#ececec', subplotpars=fig.SubplotParams(wspace=0.3, hspace=0.6, left=0.05))
	axes = set_plots(figure)
	generations_to_simulate = set_generations(generations,population_size)
	neutral_frequencies, nearly_neutral_frequencies, selection_coefficients,fitness_upper,fitness_lower = prepare_plot_data(generations_to_simulate,population_size, number_of_loci,mutation_fitness)	
	neutral_generations_to_fixation,neutral_generations_to_loss,nearly_neutral_generations_to_fixation,nearly_neutral_generations_to_loss = get_distribution_times(number_of_loci,generations_to_simulate, neutral_frequencies, nearly_neutral_frequencies)
	if not clear:
		colours = ['k','m','r','g','b']
		graph_allele_frequencies(neutral_frequencies,generations,interval, axes['ax0'],colours)
		graph_allele_frequencies(nearly_neutral_frequencies,generations,interval, axes['ax1'],colours)
		hist_time_to_fixation(neutral_generations_to_fixation,'Number of mutations','Generations until fixation',axes['ax2'])
		hist_time_to_fixation(neutral_generations_to_loss,'Number of mutations','Generations until loss',axes['ax3'])
		hist_time_to_fixation(nearly_neutral_generations_to_fixation,'Number of mutations','Generations until fixation',axes['ax4'])
		hist_time_to_fixation(nearly_neutral_generations_to_loss,'Number of mutations','Generations until loss',axes['ax5'])
		hist_selection_coefficients(np.where(nearly_neutral_generations_to_fixation > 0), selection_coefficients,fitness_lower,fitness_upper,axes['ax6'],'Number of mutations','Relative fitness \nof mutations fixed')
		hist_selection_coefficients(np.where(nearly_neutral_generations_to_loss > 0), selection_coefficients,fitness_lower,fitness_upper,axes['ax7'],'Number of mutations','Relative fitness \nof mutations lost')
	#figure.tight_layout()
	#pt.tight_layout(pad=1.5, h_pad=1.5, w_pad=1.6)
	return FigureCanvas(figure)
	
def graph_allele_frequencies(frequencies, generations,interval, axis, colours):
	plot_lines = int(generations / interval-1) #number of allele freqs to plot on allele freq by time axes
	set_to_plot = np.zeros((plot_lines,generations)) #allocate space for the random set of allele frequency trajectories to plot
	set_to_plot[0,:] = frequencies[0,0:generations]
	start_position= interval
	
	length_of_colours = len(colours)
	for i in xrange(1,plot_lines):
		for j in xrange(start_position-1, generations): #validate
			set_to_plot[i,j] = frequencies[i,j-(start_position-1)]
		start_position = start_position+interval
	
	for m in xrange(0,plot_lines):
		colour = colours[(m-1)%length_of_colours]
		axis.plot(set_to_plot[m,:],colour)
	
	axis.set_xlim(1,generations+1)
	axis.set_ylim(0.0,1.05)
	axis.set_ylabel('Allele frequency',fontdict={'size':14,'color':'#002147','family':'serif'})
	axis.set_xlabel('Generation',fontdict={'size':14,'color':'#002147','family':'serif'})	
	
	return axis
	
def hist_time_to_fixation(fixation_matrix, y_label,x_label,axis):
	mask = np.where(fixation_matrix > 0)
	axis.hist(fixation_matrix[mask])
	axis.set_ylabel(y_label,fontdict={'size':13,'color':'#002147','family':'serif'})
	axis.set_xlabel(x_label,fontdict={'size':13,'color':'#002147','family':'serif'})
	
def hist_selection_coefficients(mask, selection_coefficients,x_lim_lower,x_lim_upper,axis, y_label,x_label):
	axis.hist(selection_coefficients[mask])
	axis.set_xlim(x_lim_lower,x_lim_upper)
	axis.set_ylabel(y_label,fontdict={'size':13,'color':'#002147','family':'serif'})
	axis.set_xlabel(x_label,fontdict={'size':13,'color':'#002147','family':'serif'})
			
def set_plots(figure):
	#one column two rows
	gs0 = gridspec.GridSpec(1, 2)
	
	#four columns two rows for row 1 column 1
	gs00 = gridspec.GridSpecFromSubplotSpec(4, 4, subplot_spec=gs0[0], wspace=1.5)
	ax0 = figure.add_subplot(gs00[:-2, :])
	ax0.set_title('Strictly Neutral', fontdict={'size':18,'color':'#003E7E','family':'serif'},horizontalalignment='center', weight='bold')

	ax1 = figure.add_subplot(gs00[2:, :])	
	ax1.set_title('Nearly Neutral', fontdict={'size':18,'color':'#003E7E','family':'serif'},horizontalalignment='center', weight='bold',linespacing=3)

	#four columns two rows for row 1 column 2
	gs01 = gridspec.GridSpecFromSubplotSpec(4, 4, subplot_spec=gs0[1], wspace=1.5)
	ax2 = figure.add_subplot(gs01[0, :-2])
	ax3 = figure.add_subplot(gs01[0, -2:])
	ax4 = figure.add_subplot(gs01[-2, :-2])
	ax5 = figure.add_subplot(gs01[-2, -2:])
	ax6 = figure.add_subplot(gs01[-1, :-2])
	#ax6.xaxis.labelpad = 18
	ax7 = figure.add_subplot(gs01[-1, -2:])
	#ax7.xaxis.labelpad = 18
	
	return {'ax0':ax0,'ax1':ax1,'ax2':ax2,'ax3':ax3,'ax4':ax4,'ax5':ax5,'ax6':ax6,'ax7':ax7}
		
def prepare_plot_data(generations_to_simulate,population_size, number_of_loci, mutation_fitness):
	
	selection_coefficients, fitness_upper,fitness_lower = get_selection_coefficients(number_of_loci, mutation_fitness)
	
	neutral_frequencies = np.zeros((number_of_loci,generations_to_simulate))    # allocate space for strictly neutral allele freqs
	nearly_neutral_frequencies = neutral_frequencies.copy()    # allocate space for nearly neutral allele freqs

	neutral_frequencies[:,0] = 1/(2*population_size)  #initial allele freq at all loci since mutation is one copy
	
	nearly_neutral_frequencies[:,0] = neutral_frequencies[:,0].copy()
	
	
	for i in xrange(1,generations_to_simulate):
		neutral_frequencies = drift(neutral_frequencies, number_of_loci,population_size,i)
		nearly_neutral_frequencies = drift(nearly_neutral_frequencies,number_of_loci,population_size,i)
		nearly_neutral_frequencies = selection(nearly_neutral_frequencies,number_of_loci,i,selection_coefficients)
		
	
	
	return neutral_frequencies, nearly_neutral_frequencies, selection_coefficients, fitness_upper,fitness_lower
	

def drift(p_matrix, loci, n, gen):
	drift_matrix = p_matrix.copy()
	for k in xrange(loci): #count from 0 not 1
		if (p_matrix[k,gen-1] < 1.0) & (p_matrix[k,gen-1] > 0.0): 
			delta = (2*rand.random()-1)*sqrt((3*p_matrix[k,gen-1]*(1-p_matrix[k,gen-1]))/(2*n))
			if p_matrix[k,gen-1]+delta > 1.0:
				drift_matrix[k,gen] =1.0
			elif p_matrix[k,gen-1]+delta < 0.0:
				drift_matrix[k,gen] =0.0
			else:
				drift_matrix[k,gen] = p_matrix[k,gen-1] +delta
		else:
			drift_matrix[k,gen] =p_matrix[k,gen-1] #allele is no longer segregating so freq is constant through time
	return drift_matrix	 




	
def selection(frequencies, number_of_loci,generation, selection_coefficients):
	'''
	Allele freqs experience natural selection. Selxn model is for a
	diploid locus with additive fitness values. The selxn_coeffs fitness
	value is assumed to be for two copies of the allele, i.e. in
	homozygotes. Assumes w-bar is one. 
	Allele freq has already experienced drift. Use current generation
    allele freqs to determine change by selection.
	'''
	selection_p_matrix = frequencies.copy()
	
	for k in xrange(number_of_loci):
		if(frequencies[k,generation] < 1.0) and (frequencies[k,generation] > 0.0):
			p = frequencies[k,generation]
			q = 1-frequencies[k,generation]
			wAa = 1-0.5*(1-selection_coefficients[k,0])
			waa = selection_coefficients[k,0]
			delta = p*q*(q*(waa-wAa)+p*(wAa-1)) # get change in allele frequency due to selection
			
			if frequencies[k,generation] + delta > 1.0:
				selection_p_matrix[k,generation] = 1.0
			elif frequencies[k,generation] + delta < 0.0:
				selection_p_matrix[k,generation] = 0.0
			else:
				selection_p_matrix[k,generation] = frequencies[k,generation]+delta
		else:
			selection_p_matrix[k,generation] = frequencies[k,generation]
	return selection_p_matrix 