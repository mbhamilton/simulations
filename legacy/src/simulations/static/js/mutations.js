/* JQuery script file to be used for surveytool javascript */

	
$(document).ready(function(){
    model = $('#selectedModel').val();
    selectModel(model);
    $('#irr_mu').focus(function(){
    	selectModel(0); 
    
    });
   
    $('#irr_gens').focus(function(){
    	selectModel(0); 
    
    });
    
    $('#twf_mu').focus(function(){
    	selectModel(1); 
    
    });
    
    $('#twr_mu').focus(function(){
    	selectModel(1); 
    
    });
    
    $('#tw_gens').focus(function(){
    	selectModel(1); 
    
    });
    
    $('#nm_ne').focus(function(){
    	selectModel(2); 
    
    });
    
    $('#nm_gen').focus(function(){
    	selectModel(2); 
    
    });
    
    $('#nm_total').focus(function(){
    	selectModel(0); 
    
    });
   
function selectModel(model)
{
	if(!$('#model_'+model).attr("checked"))
		$('#model_'+model).attr("checked",true);
}
});
    
    