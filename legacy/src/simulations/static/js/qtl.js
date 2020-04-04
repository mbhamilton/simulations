/**
 * Javascript file for qtl simulation
 * Author: Marie Kolawole
 */
$(document).ready(function(){
	
	mutation_rate = $('#mr_value').text()*1000;
	initial_allele = $('#iaf_value').text()*10;
	variation_environment = $("#ve_value").text()*10;
	selection_threshold = $("#selxn_thresh_value").text()*10;
	//alert('mutation rate = '+$('#mr_value').text())
	//alert(' raw mutation rate = '+mutation_rate)
	$('#mr').val(mutation_rate);
	$('#iaf').val(initial_allele);
	$('#ve').val(variation_environment);
	$('#selxn_thresh').val(selection_threshold);
	qtl = $('#num_qtl').val();
	qtl_arch = $('#qtl_architecture').val();
	$('#qtl_'+qtl).attr("checked",true);
	$('#qtl_arch_'+qtl_arch).attr("checked",true);
	$('#nlafp_max').text(qtl)
	
	if(qtl_arch == 2)
		show_rows();
	else	
		hide_rows();
	
	$('#mr').change(function(){
		value = $(this).val();
		changeLabel('mr_value',value/1000)
	});
	$('#iaf').change(function(){
		value = $(this).val();
		changeLabel('iaf_value',value/10)
	});
	$('#ve').change(function(){
		value = $(this).val();
		changeLabel('ve_value',value/10)
	});
	$('#selxn_thresh').change(function(){
		value = $(this).val();
		changeLabel('selxn_thresh_value',value/10)
	});
	$('#max_geno').change(function(){
		value = $(this).val();
		changeLabel('max_geno_value',value)
	});
	$('#nlafp').change(function(){
		value = $(this).val();
		changeLabel('nlafp_value',value)
	});
	
	$('input[name=qtl_arch]').click(function(){
		if($(this).val()==2)
			show_rows();
		else
			hide_rows();
	});
	
	$('input[name=qtl]').click(function(){
		value = $(this).val()
		$('#nlafp').attr("max",value);
		$('#nlafp_max').text(value)
	});

});



function hide_rows(){
	$('#major_qtl_row_1').hide();
	$('#major_qtl_row_2').hide();
}

function show_rows(){
	$('#major_qtl_row_1').show();
	$('#major_qtl_row_2').show();
}