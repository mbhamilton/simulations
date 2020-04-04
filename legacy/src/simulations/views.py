# Create your views here.
from __future__ import division #using python 3.0 behaviour of division which will always return a float instead of casting to an int
from django.shortcuts import render_to_response 
from django.http import HttpResponse 
import pyximport
import numpy as np
pyximport.install(setup_args={'include_dirs':[np.get_include()]})
import _nuc_substitution as nuc
import mutation as mt
import definetti as df
import finite_H	 as fh
import qtl as qt
import mullers_ratchet as mr
import one_locus_selection as ns
import drift_selection_mutation as ds
from django.template import RequestContext
import neutral_theory as nt
#import nuc_substitution as nuc

def new_index(request):
	return render_to_response('definetti/styledindex.html',context_instance=RequestContext(request) )

def index(request):
	return render_to_response('index.html',context_instance=RequestContext(request) )
	
def definetti(request):
	model = {
			 'generationsList': xrange(0,101,10), 
			 'AAwhole':25, 
			 'Aawhole':50, 
			 'frequencysum':1-0.25-0.5,
			 'matingModel':0,
			 'fAA': 0.25,
			 'fAa': 0.5,
			 'generations':0,
			 'axis': 'checked',
			 'clear': 'true'
		 }

	updateSessionObject(request.session, model)
	return render_to_response('plot_template.html', model,context_instance=RequestContext(request) )


def clear_plot(plot_properties):
	if 'clear' in plot_properties:
		#print ' clear in session setting to true'
		clear = True
		del plot_properties['clear']
	else:
		clear = False
	return clear

def plot(request):
	#print 'plot being called\n'
	plot_properties = request.session['plot_properties']
	#print 'GET::::=== %s' %request.session['plot_properties']
	
	if 'axis' in plot_properties:
		showAxis = True
	else:
		showAxis = False
	
	# print 'prop clear before %s' %plot_properties['clear']	 
	clear = clear_plot(plot_properties)
	#print 'prop clear after %s' %plot_properties['clear']	  
	canvas=df.plot(showAxis,plot_properties['fAA'],plot_properties.get('fAa'),
				   plot_properties.get('generations'),plot_properties.get('matingModel',0),clear)
	response = HttpResponse(content_type='image/png')
	canvas.print_png(response)
	plot_properties['focusOnPlot'] = 'true'
	updateSessionObject(request.session, plot_properties)
	return response

def check_clear(request, properties):
	if 'clear' in request.GET:
		properties['clear'] = 'true'

def run(request):
	properties = request.session['plot_properties']
	properties['axis'] = request.GET.get('showAxis','off') 
	properties['fAA'] = int(request.GET.get('AA'))/100.0
	properties['fAa'] = int(request.GET.get('Aa'))/100.0
	properties['generations'] = int(request.GET.get('generationsToSimulate'))
	properties['matingModel'] = int(request.GET.get('MatingModel',0)) #in case none selected choose first
	properties['AAwhole'] = int(request.GET.get('AA'))
	properties['Aawhole'] = int(request.GET.get('Aa'))
	properties['frequencysum'] = 1-properties['fAA']-properties['fAa']
	
	if properties['axis'] =='off':
		del properties['axis'] 
	else:
		properties['axis'] ='checked'

		check_clear(request, properties)
		updateSessionObject(request.session, properties)
		return render_to_response('plot_template.html',properties,context_instance=RequestContext(request))

def updateSessionObject(session, plot_properties):
	#	check if alreay exists
	#if 'plot_properties' not in session:
	session['plot_properties']=plot_properties
	   # print 'added to session %s \n' %session['plot_properties']

def mutation(request):
	model = {
			 'irreversible_mutation': 0.0001, 
			 'irreversible_generations':100000, 
			 'twoway_forward_mutation':0.0001, 
			 'twoway_reverse_mutation':0.0005,
			 'twoway_generations':20000,
			 'population_size': 20,
			 'neutral_generations_mutation': 30,
			 'neutral_generations':200,
			 'mutation_model':0,
			 'clear':'true'
		 }
	updateSessionObject(request.session, model)
	return render_to_response('mutation_template.html',model,context_instance=RequestContext(request) )

def run_mutation(request):
	properties = request.session['plot_properties']
	mutation_model = int(request.GET.get('mutationModel',0))
	
	check_clear(request, properties)
	
	if mutation_model == 0:
		mu = mt.validate_mu(float(request.GET.get('irr_mu',0.0001)))
		properties['irreversible_mutation']=mu
		gens = mt.validate_irreversible_generations(int(request.GET.get('irr_gens',100000)))
		properties['irreversible_generations']=gens
	   # print ' arguments mu=', mu,'gens=', gens
		
	elif mutation_model == 1:
	   # print 'what is mu? %s' %request.GET.get('twf_mu')
		mu = mt.validate_mu(float(request.GET.get('twf_mu',0.0001)))
	   # print 'what is mu? %f' %mu
		properties['twoway_forward_mutation']=mu
		nu = mt.validate_mu(float(request.GET.get('twr_mu',0.0005)))
		properties['twoway_reverse_mutation']=nu
		gens = mt.validate_twoway_generations(int(request.GET.get('tw_gens',20000)))
		properties['twoway_generations']= gens
		
	elif mutation_model == 2:
		ne = mt.validate_ne(int(request.GET.get('nm_ne',20)))
		interval = mt.validate_generations_between(int(request.GET.get('nm_gen',30)))
		gens = mt.validate_total_generations(int(request.GET.get('nm_total',200)))
		properties['population_size']=ne
		properties['neutral_generations_mutation']=interval
		properties['neutral_generations']=gens
	
	properties['mutation_model']= mutation_model
	updateSessionObject(request.session, properties)
	return render_to_response('mutation_template.html',properties,context_instance=RequestContext(request))


def plot_mutation(request):
	
	properties = request.session['plot_properties']
	mutation_model = properties['mutation_model']
	clear =	 clear_plot(properties)
	if mutation_model == 0:
		canvas = mt.plot(mutation_model,clear,properties['irreversible_generations'], properties['irreversible_mutation'])
	elif mutation_model == 1:
		canvas = mt.plot(mutation_model,clear,properties['twoway_generations'], properties['twoway_forward_mutation'],properties['twoway_reverse_mutation'])
	elif mutation_model == 2:
		canvas = mt.plot_neutral_mutation(properties['population_size'],properties['neutral_generations'], properties['neutral_generations_mutation'],clear)
		
		
	response = HttpResponse(content_type='image/png')
	canvas.print_png(response)
	properties['focusOnPlot'] = 'true'
	updateSessionObject(request.session, properties)
	return response

def finiteh(request):
	model = {
			 'populationSize': 50, 
			 'alleleFrequency':5, #0.5 
			 'numberOfGenerations':100, 
			 'replicatePopulations': 6,
			 'clear': 'true'
		 }
	updateSessionObject(request.session, model)
	return render_to_response('finite_H_template.html',model,context_instance=RequestContext(request) )

def plot_finiteh(request):
	#print 'plot finite request map %s' %request.session['plot_properties']
	properties = request.session['plot_properties']
	populationSize = properties['populationSize']
	numberOfGenerations = properties['numberOfGenerations']
	alleleFrequency = properties['alleleFrequency']
	replicatePopulations = properties['replicatePopulations']
	clear = clear_plot(properties)
	canvas = fh.drift_h(populationSize, numberOfGenerations, alleleFrequency/10,replicatePopulations, clear)
	response = HttpResponse(content_type='image/png')
	canvas.print_png(response)
	properties['focusOnPlot'] = 'true'
	updateSessionObject(request.session, properties)
	return response
	
def run_finiteh(request):
	properties = request.session['plot_properties']
	properties['populationSize'] = fh.validate_population_size(int(request.GET.get('populationSize',50)))
	properties['numberOfGenerations'] = fh.validate_number_of_generations(int(request.GET.get('numberOfGenerations',100)))
	properties['alleleFrequency'] = int(request.GET.get('alleleFrequency',5))
	properties['replicatePopulations'] = int(request.GET.get('replicatePopulations',6))
	check_clear(request, properties)
	updateSessionObject(request.session, properties)
	return render_to_response('finite_H_template.html',properties,context_instance=RequestContext(request) )

def mullers(request):
	model = {
			 'populationSize': 200, 
			 'generations':200, 
			 'selection':0.01, 
			 'mutation':0.06,
			 'loci':100,
			 'mean': 1,
			 'max': 10, 
			 'clear':'true'
		 }
	updateSessionObject(request.session, model)
	return render_to_response('mullers_ratchet_template.html',model,context_instance=RequestContext(request) )

def plot_mullers(request):
	#print 'plot mullers request map %s' %request.session['plot_properties']
	properties = request.session['plot_properties']
	populationSize = properties['populationSize']
	generations = properties['generations']
	selection = properties['selection']
	mutation = properties['mutation']
	loci = properties['loci']
	mean = properties['mean']
	maximum = properties['max']
	clear = clear_plot(properties)
	canvas = mr.plot(populationSize, generations, loci,selection,mutation,mean,maximum, clear)
	response = HttpResponse(content_type='image/png')
	canvas.print_png(response)
	properties['focusOnPlot'] = 'true'
	updateSessionObject(request.session, properties)
	return response

def run_mullers(request):
	properties = request.session['plot_properties']
	properties['populationSize'] = mr.validate_population_size(int(request.GET.get('populationSize',200)))
	properties['generations'] = mr.validate_generations(int(request.GET.get('generations',200)))
	properties['selection'] =  mr.validate_coefficient(float(request.GET.get('selection',0.01)))
	properties['mutation'] =  mr.validate_mutation(float(request.GET.get('mutation',0.06)))
	properties['loci'] = mr.validate_loci(int(request.GET.get('loci',100)))
	properties['mean'] =  mr.validate_mean_family_size(float(request.GET.get('mean',1)))
	properties['max'] =	 mr.validate_max_family_size(int(request.GET.get('max',10))) 
	check_clear(request, properties)
	updateSessionObject(request.session, properties)
	return render_to_response('mullers_ratchet_template.html',properties,context_instance=RequestContext(request) )

def naturalselection(request):
	model = {
			 'pAA': 0.25, 
			 'pAa': 0.5, 
			 'paa': 0.25, 
			 'wAA': 1.2,
			 'wAa': 1.1,
			 'waa': 1.0,
			 'gens': 100, 
			 'clear':'true'
			 
		 }
	updateSessionObject(request.session, model)
	return render_to_response('natural_selection_template.html',model,context_instance=RequestContext(request))

def plot_naturalselection(request):
	properties = request.session['plot_properties']
	pAA = properties['pAA']
	pAa = properties['pAa']
	paa = properties['paa']
	wAA = properties['wAA']
	wAa = properties['wAa']
	waa = properties['waa']
	generations = properties['gens']
	clear = clear_plot(properties)
	canvas = ns.plot(pAA,pAa,paa,wAA,wAa,waa,generations, clear)
	response = HttpResponse(content_type='image/png')
	canvas.print_png(response)
	properties['focusOnPlot'] = 'true'
	updateSessionObject(request.session, properties)
	return response

def run_naturalselection(request):
	properties = request.session['plot_properties']
	
	properties['pAA'] = float(request.GET.get('pAA',25))/100
	properties['pAa'] = float(request.GET.get('pAa',50))/100
	properties['paa'] = 1- properties['pAA']-properties['pAa']
	properties['wAA'] = float(request.GET.get('wAA',12))/10
	properties['wAa'] = float(request.GET.get('wAa',11))/10
	properties['waa'] = float(request.GET.get('waa',10))/10
	properties['gens'] =  ns.validate_generations(int(request.GET.get('gens',100))) 
	check_clear(request, properties)
	updateSessionObject(request.session, properties)
	return render_to_response('natural_selection_template.html',properties,context_instance=RequestContext(request) )

def qtl(request):
	model = {
			 'ne': 20, 
			 'mutation_rate': 0.001, 
			 'initial_allele': 0.1, 
			 'variation_environment': 0.1,
			 'max_genotype': 10,
			 'selxn_thresh': 0.2,
			 'generations': 100, 
			 'number_of_major_qtl': 2,
			 'number_of_allele_frequenccy': 10,
			 'genotypic_variation':0.2,
			 'qtl':10,
			 'nplotloci': 10,
			 'qtl_arch': 1,
			 'clear':'true'
			 
		 }
	updateSessionObject(request.session, model)
	return render_to_response('qtl_template.html',model,context_instance=RequestContext(request))

def plot_qtl(request):
	properties = request.session['plot_properties']
	populationSize = properties['ne']
	mu = properties['mutation_rate']
	initial_allele = properties['initial_allele']
	variation_environment = properties['variation_environment']
	max_genotype = properties['max_genotype']
	selxn_thresh = properties['selxn_thresh']
	number_of_major_qtl = properties['number_of_major_qtl']
	generations = properties['generations']
	genotypic_variation = properties['genotypic_variation']
	number_of_qtl = properties['qtl']
	nplotloci = properties['nplotloci']
	qtl_arch = properties['qtl_arch']
	
	clear = clear_plot(properties)
	canvas = qt.plot(populationSize,number_of_qtl,mu,generations,initial_allele,variation_environment,max_genotype, selxn_thresh,qtl_arch,genotypic_variation,number_of_major_qtl,nplotloci, clear)
	response = HttpResponse(content_type='image/png')
	canvas.print_png(response)
	properties['focusOnPlot'] = 'true'
	updateSessionObject(request.session, properties)
	return response

def run_qtl(request):
	properties = request.session['plot_properties']
	properties['ne'] = qt.validate_population_size(int(request.GET.get('ne',20)))
	properties['mutation_rate'] = float(request.GET.get('mr',1))/1000
	properties['initial_allele'] = float(request.GET.get('iaf',1))/10
	properties['variation_environment'] = float(request.GET.get('ve',1))/10
	properties['max_genotype'] = int(request.GET.get('max_geno',10))
	properties['selxn_thresh'] = float(request.GET.get('selxn_thresh',2))/10
	properties['number_of_major_qtl'] = qt.validate_number_of_major_qtl(int(request.GET.get('num_major_qtl',2)))
	properties['generations'] = qt.validate_generations(int(request.GET.get('generations',100)))
	properties['genotypic_variation'] = qt.validate_genotypic_variance(float(request.GET.get('geno_var',0.2)),properties['number_of_major_qtl'])
	properties['qtl'] = int(request.GET.get('qtl',10))
	properties['nplotloci'] = int(request.GET.get('nlafp',10))
	properties['qtl_arch'] = int(request.GET.get('qtl_arch',1))
	
	check_clear(request, properties)
	updateSessionObject(request.session, properties)
	return render_to_response('qtl_template.html',properties,context_instance=RequestContext(request) )

def driftselection(request):
	model = {
			 'replicate_population': 200,
			 'ne': 20, 
			 'fa': 0.5,
			 'number_of_generations': 100, 
			 'rpia':10,
			 'wAA':1.2,
			 'wAa':1.1,
			 'waa':1.0,
			 'mu':0.01,
			 'nu':0.01,
			 'natural_selection':'on',
			 'clear':'true'
			 
		 }
	updateSessionObject(request.session, model)
	return render_to_response('drift_selection_template.html',model,context_instance=RequestContext(request))

def plot_driftselection(request):
	properties = request.session['plot_properties']
	number_of_replicate_populations = properties['replicate_population']
	populationSize = properties['ne']
	initial_allele = properties['fa']
	number_of_generations = properties['number_of_generations']
	number_of_replicate_populations_in_allele = properties['rpia']
	wAA = properties['wAA']
	wAa = properties['wAa']
	waa = properties['waa']
	mu = properties['mu']
	nu = properties['nu']
	
	if 'natural_selection' in properties:
		natural_selection = True
	else:
		natural_selection = False
	
	if 'mutation' in properties:
		mutation = True
	else:
		mutation = False
	clear = clear_plot(properties)
	canvas = ds.plot(number_of_replicate_populations,populationSize,initial_allele,number_of_generations, number_of_replicate_populations_in_allele,natural_selection,wAA,wAa,waa,mutation,mu,nu, clear)
	response = HttpResponse(content_type='image/png')
	canvas.print_png(response)
	properties['focusOnPlot'] = 'true'
	updateSessionObject(request.session, properties)
	return response

def run_driftselection(request):
	properties = request.session['plot_properties']
	
	properties['replicate_population'] = ds.validate_number_of_replicate_populations(int(request.GET.get('replicate_population',200)))
	properties['ne'] = ds.validate_value(int(request.GET.get('ne',20)))
	properties['fa'] = float(request.GET.get('initial_allele_frequency',5))/10
	properties['number_of_generations'] = int(request.GET.get('number_of_generations',100))
	properties['rpia'] = int(request.GET.get('replicate_populations_in _allele',10))
	properties['natural_selection']= request.GET.get('natural_selection','off')
	
	if properties['natural_selection'] =='off':
		del properties['natural_selection'] 
	else:
		properties['natural_selection'] ='checked'
		
	properties['wAA'] = float(request.GET.get('wAA',12))/10
	properties['wAa'] = float(request.GET.get('wAa',11))/10
	properties['waa'] = float(request.GET.get('waa',10))/10
	properties['mutation'] = request.GET.get('mutation','off')
	if properties['mutation'] =='off':
		del properties['mutation'] 
	else:
		properties['mutation'] ='checked'

	properties['mu'] =	float(request.GET.get('mu',10))/1000 
	properties['nu'] =	float(request.GET.get('nu',10))/1000 
	
	check_clear(request, properties)
	updateSessionObject(request.session, properties)
	return render_to_response('drift_selection_template.html',properties,context_instance=RequestContext(request))
	
def neutral_theory(request):
	model = {
			'generations_to_view':200,
			'effective_population':20,
			'generations_between_mutations':20,
			'number_of_mutations': 500,
			'relative_fitness':5,
			'clear':'true'	
	}
	updateSessionObject(request.session, model)
	return render_to_response('neutral_theory_template.html',model,context_instance=RequestContext(request))
	
def run_neutral_theory(request):
	properties = request.session['plot_properties']
	properties['generations_to_view'] = nt.validate_generations(int(request.GET.get('gensToView',200)))
	properties['effective_population'] = nt.validate_population_size(int(request.GET.get('populationSize',20)))
	properties['generations_between_mutations'] = nt.validate_generations_between_mutations(int(request.GET.get('gensBetweenMuts',20)))
	properties['number_of_mutations'] = nt.validate_number_of_mutations(int(request.GET.get('numberOfLoci',500)))
	properties['relative_fitness'] = int(request.GET.get('relativeFitness',50))/10
	check_clear(request, properties)
	updateSessionObject(request.session, properties)
	return render_to_response('neutral_theory_template.html',properties,context_instance=RequestContext(request))
	
def plot_neutral_theory(request):
	properties = request.session['plot_properties']
	clear = clear_plot(properties)
	canvas = nt.plot(properties['generations_to_view'],properties['effective_population'],properties['number_of_mutations'] , properties['relative_fitness'],properties['generations_between_mutations'], clear)
	response = HttpResponse(content_type='image/png')
	canvas.print_png(response)
	properties['focusOnPlot'] = 'true'
	updateSessionObject(request.session, properties)
	return response

def nuc_substitution(request):
	model = {
			'a':0.25,
			't':0.25,
			'c':0.25,
			'g':0.25,
			'bp':2500,
			'gens':500,
			'ag':0.0005,
			'ct':0.0005,
			'ac':0.0005,
			'at':0.0005,
			'gt':0.0005,
			'gc':0.0005,
			'clear':'true'	
	}
	updateSessionObject(request.session, model)
	return render_to_response('nuc_substitution_template.html',model,context_instance=RequestContext(request))

def run_nuc_substitution(request):
	properties = request.session['plot_properties']
	properties['ag'] = nuc.validate_transition(float(request.GET.get('ag',0.0005)))
	properties['ct'] = nuc.validate_transition(float(request.GET.get('ct',0.0005)))
	properties['ac'] = nuc.validate_transition(float(request.GET.get('ac',0.0005)))
	properties['at'] = nuc.validate_transition(float(request.GET.get('at',0.0005)))
	properties['gt'] = nuc.validate_transition(float(request.GET.get('gt',0.0005)))
	properties['gc'] = nuc.validate_transition(float(request.GET.get('gc',0.0005)))
	
	properties['a'] = float(request.GET.get('A',0.25))
	properties['t'] = float(request.GET.get('T',0.25))
	properties['c'] = float(request.GET.get('C',0.25))
	properties['g'] = nuc.get_G(properties['a'], properties['c'], properties['t'])
	properties['bp'] = int(request.GET.get('dna',2500))
	properties['gens'] = int(request.GET.get('generations',500))
	properties['exampleModel'] = request.GET.get('exampleModel','')
	check_clear(request, properties)
	updateSessionObject(request.session, properties)
	return render_to_response('nuc_substitution_template.html',properties,context_instance=RequestContext(request))

def plot_nuc_substitution(request):
	properties = request.session['plot_properties']
	clear = clear_plot(properties)
	canvas = nuc.plot(properties['a'], properties['t'], properties['c'], properties['g'], properties['ag'], properties['ct'], properties['ac'], properties['at'], properties['gc'], properties['gt'], properties['bp'], properties['gens'], clear)
	response = HttpResponse(content_type='image/png')
	canvas.print_png(response)
	properties['focusOnPlot'] = 'true'
	updateSessionObject(request.session, properties)
	return response
