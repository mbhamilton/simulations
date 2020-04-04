/* JQuery script file to be used for surveytool javascript */

	
$(document).ready(function(){
	divisor=10;

    $('#replicatePopulationsValue').text(
    		$('#replicatePopulations').val()
    		);
    $('#alleleFrequencyValue').text(
    		$('#alleleFrequency').val()/divisor);
    
    $('#alleleFrequency').change(function(){
    	$('#alleleFrequencyValue').text(
        		$('#alleleFrequency').val()/divisor);
    });
    
    $('#replicatePopulations').change(function(){
    	$('#replicatePopulationsValue').text(
        		$('#replicatePopulations').val());
    });
});
    
    