/* JQuery script file to be used for surveytool javascript */

$(document).ready(function(){
   //alert('checking generations '+$('#generations'))
    model = $('#selectedModel').val()
    generation = $('#selectedgenerations').val()
    $('#MatingModel_'+model).attr("checked",true);
    
    $('select[name=generationsToSimulate] option[value='+generation+']').attr('selected',true);
    
    $('#AA').change(function(){
        fAA =parseInt($(this).val())
        fAa = parseInt($('#Aa').val())
        fAA = ensure_geneotype_frequency_sum_is_one('AA',fAa,fAA)
        changeValue('fAA','aa',fAA/100,(100-fAA-fAa)/100);
        
        
        });
    $('#Aa').change(function(){
        fAA = parseInt($('#AA').val())
        fAa = parseInt($(this).val())
        fAa = ensure_geneotype_frequency_sum_is_one('Aa',fAa,fAA)
        changeValue('fAa','aa',fAa/100,(100-fAA-fAa)/100);
        });
    
$("form").submit(function() {

    $('#run').text('Processing...')
    });

    });
