/**
 * Javascript file for qtl simulation
 * Author: Marie Kolawole
 */
$(document).ready(function(){
	
	fitness = $('#fit_value').text()*10;
//	alert(fitness)
	$('#relfit').val(fitness);
//	alert($('#relfit').val())
	
	$('#relfit').change(function(){
		value = $(this).val();
		changeLabel('fit_value',value/10)
	});

	
});

