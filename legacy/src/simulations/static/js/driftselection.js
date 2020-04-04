/**
 * Javascript file for drift selection  simulation
 * uses simulationFunctions file function
 * Author: Marie Kolawole
 */
$(document).ready(function(){
	
	fa = $('#fa_value').text()*10;
	rpia = $('#rpia_value').text();
	wAA = $("#wAA_value").text()*10;
	wAa = $("#wAa_value").text()*10;
	waa = $("#waa_value").text()*10;
	mu = $("#mu_value").text()*1000;
	nu = $("#nu_value").text()*1000;
	//alert('mutation rate = '+$('#mr_value').text())
	//alert(' raw mutation rate = '+mutation_rate)
	$('#initial_allele_frequency').val(fa);
	$('#replicate_populations_in_allele').val(rpia);
	$('#wAA').val(wAA);
	$('#wAa').val(wAa);
	$('#waa').val(waa);
	$('#mu').val(mu);
	$('#nu').val(nu);
	

	$('#initial_allele_frequency').change(function(){
		value = $(this).val();
		changeLabel('fa_value',value/10)
	});
	$('#replicate_populations_in_allele').change(function(){
		value = $(this).val();
		changeLabel('rpia_value',value)
	});
	$('#wAA').change(function(){
		value = $(this).val();
		changeLabel('wAA_value',value/10)
	});
	$('#wAa').change(function(){
		value = $(this).val();
		changeLabel('wAa_value',value/10)
	});
	$('#waa').change(function(){
		value = $(this).val();
		changeLabel('waa_value',value/10)
	});
	$('#mu').change(function(){
		value = $(this).val();
		changeLabel('mu_value',value/1000)
	});
	$('#nu').change(function(){
		value = $(this).val();
		changeLabel('nu_value',value/1000)
	});
	
});
