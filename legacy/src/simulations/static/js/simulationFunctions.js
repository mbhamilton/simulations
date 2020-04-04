/**
 * Common functions used in javascript files for simulation
 */



//Functions
function changeValue(sliderlabel, aalabel, value, sum){
    $('#'+sliderlabel).text(Math.round(value*100)/100); //2 decimal places
// rounding to two decimal places, some instances i.e 0.56 of Aa where giving about 17 decimal places
    $('#'+aalabel).text(Math.round(sum*100)/100); 
}

function changeLabel(sliderlabel, value){
	 $('#'+sliderlabel).text(value);//(Math.round(value*100)/100); //2 decimal places
}
 
function ensure_geneotype_frequency_sum_is_one(sliderid,fAa,fAA){
    // alert('fAA '+fAA+" , fAa "+fAa)
     //ensure that the sum of the genotype frequencies is one or less
     //write unit test for this method
 	
     frequencyValue=0
     if (sliderid.match('Aa$')){
    	 //alert('match Aa')
         frequencyValue = fAa/100
         max_fAa = 1- (fAA/100)
         if (frequencyValue > max_fAa)
             frequencyValue = max_fAa 
     }
     else if (sliderid.match('AA$')){
    	 //alert('match AA')
         frequencyValue= fAA/100
         max_fAA = 1- (fAa/100)
         if (frequencyValue > max_fAA)
             frequencyValue = max_fAA
     }
     frequencyValue*=100 //make whole number
     //alert('#'+sliderid)
     $('#'+sliderid).val(frequencyValue)
     return frequencyValue
 }

//shorthand for $(document).ready()
$(function(){
	if($('#focusOnPlot').val()){
		if($('#fullplot').val() != undefined){
			element = $('#fullplot')
		}
		else if($('#plot').val() != undefined){
			element= $('#plot');
		}
				
		element.css('outline','none').attr('tabindex',-1).focus();
	
	}

});