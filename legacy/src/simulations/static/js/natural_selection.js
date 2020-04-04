/**
 * uses functions in simulationFunctions.js
 */$(document).ready(function(){
   //alert('checking generations '+$('#generations'))
    
	 pAA = $('#pAA').text()*100
	 pAa = $('#pAa').text()*100
	 //alert(pAa)
	 $('#AA').val(pAA)
	 $('#Aa').val(pAa)
	 wAA = $('#wAA_value').text()*10
	 wAa = $('#wAa_value').text()*10
	 waa = $('#waa_value').text()*10
	 $('#wAA').val(wAA)
	 $('#wAa').val(wAa)
	 $('#waa').val(waa)
    
	  $('#AA').change(function(){
        pAA =parseInt($(this).val())
        pAa = parseInt($('#Aa').val())
        pAA = ensure_geneotype_frequency_sum_is_one('AA',pAa,pAA)
        changeValue('pAA','paa',pAA/100,(100-pAA-pAa)/100);
        });
	 
    $('#Aa').change(function(){
        pAA = parseInt($('#AA').val())
        pAa = parseInt($(this).val())
        pAa = ensure_geneotype_frequency_sum_is_one('Aa',pAa,pAA)
        changeValue('pAa','paa',pAa/100,(100-pAA-pAa)/100);
        });
    
    $('#wAA').change(function(){
        wAA =parseInt($(this).val())
        $('#wAA_value').text(wAA/10)        
        });
    
    $('#wAa').change(function(){
    	 wAa =parseInt($(this).val())
         $('#wAa_value').text(wAa/10)
         });
    
    $('#waa').change(function(){
    	 waa =parseInt($(this).val())
         $('#waa_value').text(waa/10)
         });

    });
