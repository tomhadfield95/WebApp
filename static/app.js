let coordStore = []
let filesToUpload = []
let filesToUploadFull = []

function uploadFile(form){ 
	const formData = new FormData(form);
	let oOutput = document.getElementById("static_file_response")
	var oReq = new XMLHttpRequest();
	oReq.open("POST", "upload_static_file", true);
	oReq.onload = function(oEvent) {
     		if (oReq.status == 200) {
       			oOutput.innerHTML = "Successfully uploaded ligand file!";
       			console.log(oReq.response)
     		} else {
       			oOutput.innerHTML = "Error occurred when trying to upload your file.<br \/>";
     		}
     	}; 

	oOutput.innerHTML = "Sending file!";
 	console.log("Sending file!")
	oReq.send(formData);

}




function uploadMultipleFiles(form){ 
	const formData = new FormData(form);
	let oOutput = document.getElementById("static_file_response")
	var oReq = new XMLHttpRequest();
	oReq.open("POST", "upload_multiple_static_files", true);
	oReq.onload = function(oEvent) {
     		if (oReq.status == 200) {
       			oOutput.innerHTML = "Successfully uploaded files";
       			console.log(oReq.response)
       			//console.log(JSON.parse(oReq.response).files_to_upload)
       			filesToUpload = filesToUpload.concat(JSON.parse(oReq.response).files_to_upload)
       			console.log('Files to upload variable')
       			console.log(filesToUpload)
     		} else {
       			oOutput.innerHTML = "Error occurred when trying to upload your file.<br \/>";
     		}
     	}; 

	oOutput.innerHTML = "Sending file!";
 	console.log("Sending file!")
	oReq.send(formData);

}




function displayMol(img_id, input_id){
    console.log(img_id)
    console.log(document.body.getElementsByTagName("*"))
    console.log(document.getElementById(input_id).files[0].name)
    console.log("/static/" + document.getElementById(input_id).files[0].name)
    //document.getElementById(img_id).src = "/static/imgs/" + document.getElementById(input_id).files[0].name
    document.getElementById(img_id).src = "/static/imgs/mol_with_indices.png"
    
    $("#select_index").removeClass('hidden')
}


function displayLinkedMols(webapp_session_directory, img_id){
    //console.log(img_id)
    //console.log(document.body.getElementsByTagName("*"))
    //console.log(document.getElementById(input_id).files[0].name)
    //console.log("/static/" + document.getElementById(input_id).files[0].name)
    //document.getElementById(img_id).src = "/static/imgs/" + document.getElementById(input_id).files[0].name
    //document.getElementById(img_id).src = "/static/imgs/linked_mols.svg"
    document.getElementById(img_id).src = webapp_session_directory + '/linked_mols.svg'
    
    //$("#select_index").removeClass('hidden')
}


function displayDeLinkerInput(webapp_session_directory, img_id){
    
    document.getElementById("check_input_div").style.display="block"
    document.getElementById(img_id).src = webapp_session_directory + '/DeLinker_input.svg'
    
 
}




function display3DMol(){

        coordStore.length = 0 //delete any existing coordinates which have been saved
        
        let element = $('#mol_3d');
        let config = { backgroundColor: 'white' };
        let viewer = $3Dmol.createViewer( element, config );
        let pdbUri = '/static/imgs/ligand.sdf';
        //let pdbUri="{{ upload_folder }} + '/ligand1.sdf'"
        
       
        console.log(pdbUri)
        jQuery.ajax( pdbUri, { 
            success: function(data) {
            
              $("#instructions").hide();
              let v = viewer;
              v.addModel( data, "sdf" );                       /* load data */
              v.setStyle({stick:{radius:0.15}});               /* style all atoms */
              
              v.setClickable({},true,function(atom,viewer,event,container) { 
                        v.addSphere({center: {x:atom.x, y:atom.y, z:atom.z}, radius: 0.5, color:'red'})
                        v.render()
                        coordStore.push([atom.x, atom.y, atom.z])     
                   });
                   
              v.zoomTo();                                      /* set camera */
              v.render();                                      /* render scene */
              v.zoom(1.2, 1000);                               /* slight zoom */
            },
            error: function(hdr, status, err) {
              console.error( "Failed to load PDB " + pdbUri + ": " + err );
            },
        });
    
    }




function display3DMol_2(molDir){

        coordStore.length = 0 //delete any existing coordinates which have been saved
        
        let element = $('#mol_3d');
        let config = { backgroundColor: 'white' };
        let viewer = $3Dmol.createViewer( element, config );
        //let pdbUri = '/static/imgs/ligand.sdf';
        //let pdbUri="{{ upload_folder }} + '/ligand1.sdf'"
        
        
        //Get full file path of molecules we want to view
        for (let i = 0; i < filesToUpload.length; i++){
            filesToUploadFull.push(molDir + "/" +filesToUpload[i])
        }
        
        console.log(filesToUploadFull)
        
        
        for (let i = 0; i < filesToUploadFull.length; i++){

            let pdbUri = filesToUploadFull[i]        
            console.log(pdbUri)
            
            if (i == 0){
                //Get rid of instructions
                
                $("#instructions").hide();
                document.getElementById("user_params_form").style.display="block";
            }
            
            
            if (pdbUri.includes("ligand")){
            
            jQuery.ajax( pdbUri, { 
            success: function(data) {
                viewer.addModel( data, "sdf" );                       /* load data */
                viewer.setStyle({stick:{radius:0.15}});               /* style all atoms */    
                   
                viewer.setClickable({},true,function(atom,viewer,event,container) { 
                    console.log('Clicked')
                    viewer.addSphere({center: {x:atom.x, y:atom.y, z:atom.z}, radius: 0.5, color:'red'})
                    viewer.render()
                    coordStore.push([atom.x, atom.y, atom.z])     
                   });
                   
                viewer.zoomTo();                                      /* set camera */
                viewer.render();                                      /* render scene */
                viewer.zoom(1.2, 1000);                               /* slight zoom */
                
                
                console.log('Here')
                
                
                },
            error: function(hdr, status, err) {
              console.error( "Failed to load file " + pdbUri + ": " + err );
                },
            });
            
                
            } else {
                
            jQuery.ajax( pdbUri, { 
            success: function(data) {
                viewer.addModel( data, "pdb" );                       /* load data */
                viewer.zoomTo();                                      /* set camera */
                viewer.render();                                      /* render scene */
                viewer.zoom(1.2, 1000);                               /* slight zoom */
                },
            error: function(hdr, status, err) {
              console.error( "Failed to load file " + pdbUri + ": " + err );
                },
            });
            
            }
            
            
            

            
       
        } //end for loop
    } //end function 




function printSavedCoords(){
    console.log(coordStore)
}


$(document).ready(function () {
    $("#save_coords").on("click", function() {
        var js_data = JSON.stringify(coordStore);
        $.ajax({                        
            url: '/save_array',
            type : 'post',
            contentType: 'application/json',
            dataType : 'json',
            data : js_data
        }).done(function(result) {
            console.log(result);
            /*$("#data").html(result);*/
            $("#data").html('Successfully saved coordinates!')
        }).fail(function(jqXHR, textStatus, errorThrown) {
            console.log("fail: ",textStatus, errorThrown);
        });
    });
});


$(document).ready(function() {
    $("#generate_elabs").on("click", function() {
    
        $("#text_generating").attr("style", "display:block")
    
    
        var js_data = JSON.stringify(coordStore);
        $.ajax({                        
            url: '/generate_linkers',
            type : 'post',
            contentType: 'application/json',
            dataType : 'json',
            data : js_data
        }).done(function(result) {
           
            /*console.log(result);*/
            /*$("#data").html(result);*/
            $("#data").html('Successfully generated linkers!')
             $("#text_generating").attr("style", "display:none")
             $("#finished_generating").attr("style", "display:block")
            
        }).fail(function(jqXHR, textStatus, errorThrown) {
            console.log("fail: ",textStatus, errorThrown);
        });
    });
});





$(document).ready(function() {
  // Wait for the content to load.
  $("form[name='input_form']").submit(function(evt) {
    // If the form is to be submitted, ignore the standard behavior.
    evt.preventDefault();
    // Serialize the inputs to an array.
    let inputFields = $(this).serializeArray();
    // Send the data as JSON via AJAX.
    $.ajax({
      method: "POST",
      url: "/save_job_parameters",
      contentType: "application/json;charset=utf-8",
      dataType: "json",
      data: JSON.stringify({ input_fields: inputFields })
    }).done(data => {
      // Use the response here.
      console.log(data);
    });
  });
});
    








  /*
            
            */

