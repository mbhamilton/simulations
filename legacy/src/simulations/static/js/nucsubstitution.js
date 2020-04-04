/* JQuery script file to be used for surveytool javascript */

$(document).ready(function(){
   //alert("Nothing!!")
	//alert('Focus on plot '+$('#generations'))
    model = $('#selectedModel').val()
    $('#'+model).attr("checked",true);
        
    $("input[name='exampleModel']").change(function(){
    	var model = $(this).val()
        if( model == 'jc69'){
        	changeFrequenciesOrRates(0.25, ['#A','#T','#C','#G']);
        	changeFrequenciesOrRates(0.0005,['#ag','#ct','#ac','#at','#gc','#gt']);	
        	
        }
        else if ( model =='k80'){
        	changeFrequenciesOrRates(0.25, ['#A','#T','#C','#G']);
        	changeFrequenciesOrRates(0.0005,['#ag','#ct']);
        	changeFrequenciesOrRates(0.00001,['#ac','#at','#gc','#gt']);
        }
        else if ( model =='t92'){
        	changeFrequenciesOrRates(0.18, ['#A','#T']);
        	changeFrequenciesOrRates(0.32, ['#C','#G']);
        	changeFrequenciesOrRates(0.0005,['#ag','#ct']);
        	changeFrequenciesOrRates(0.00001,['#ac','#at','#gc','#gt']);
        }
        else if (model =='tn93'){
        	changeFrequenciesOrRates(0.18, ['#A','#T']);
        	changeFrequenciesOrRates(0.32, ['#C','#G']);
        	changeFrequenciesOrRates(0.0005,['#ag']);
        	changeFrequenciesOrRates(0.00001,['#ct','#ac','#at','#gc','#gt']);
        }
        });
    
    
    function changeFrequenciesOrRates(freq, list){
    	
    	$.each(list,function(index,value){
    		if($(value).is('input'))
    			$(value).val(freq);
    		else if($(value).is('label'))
    			$(value).text(freq);
    	});
    }
    
    $('#dna').change(function(){
		value = $(this).val();
		changeLabel('dnaValue',value);
	});
    
    $('#generations').change(function(){
		value = $(this).val();
		changeLabel('genValue',value);
	});
    
    $('#A').change(function(){
        a = parseFloat($(this).val());
        c = parseFloat($('#C').val());
        t = parseFloat($('#T').val());
        min=0.0;
        max=1-c-t;
        alert('value entered '+a)
        if(a<min)
        	a=min;
        else if(a>max)
        	a = max;
        
        g = Math.abs(1-(a+c+t))
        $(this).val(a);
        changeLabel('G',Math.round(g*100)/100);
        });
    
    $('#T').change(function(){
    	t = parseFloat($(this).val());
        c = parseFloat($('#C').val());
        a = parseFloat($('#A').val());
        min=0.0;
        max=1-c-a;
        //alert('value entered '+t)
        if(t<min)
        	t=min;
        else if(t>max)
        	t = max;
        
        g = Math.abs(1-(a+c+t))
        $(this).val(t);
        changeLabel('G',Math.round(g*100)/100);
        });
    
    $('#C').change(function(){
    	c = parseFloat($(this).val());
        a = parseFloat($('#A').val());
        t = parseFloat($('#T').val());
        min=0.0;
        max=1-t-a;
       // alert('value entered '+c)
        if(c<min)
        	c=min;
        else if(c>max)
        	c = max;
        
        g = Math.abs(1-(a+c+t))
        $(this).val(c);
        //round to 2 decimal places
        changeLabel('G',Math.round(g*100)/100);
        });
    
    $("form").submit(function() {

    $('#run').val('Processing...')
    });

    });
