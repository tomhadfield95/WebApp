let coordStore = []

function uploadFile(form){ 
	const formData = new FormData(form);
	let oOutput = document.getElementById("static_file_response")
	var oReq = new XMLHttpRequest();
	oReq.open("POST", "upload_static_file", true);
	oReq.onload = function(oEvent) {
     		if (oReq.status == 200) {
       			oOutput.innerHTML = "Successfuly uploaded ligand file!";
       			console.log(oReq.response)
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


function display3DMol(){

        coordStore.length = 0 //delete any existing coordinates which have been saved
        
        let element = $('#mol_3d');
        let config = { backgroundColor: 'white' };
        let viewer = $3Dmol.createViewer( element, config );
        let pdbUri = '/static/imgs/ligand.sdf';
        jQuery.ajax( pdbUri, { 
            success: function(data) {
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


function printSavedCoords(){
    console.log(coordStore)
}


$(document).ready(function () {
    $("#write_coords").on("click", function() {
        var js_data = JSON.stringify(coordStore);
        $.ajax({                        
            url: '/save_array',
            type : 'post',
            contentType: 'application/json',
            dataType : 'json',
            data : js_data
        }).done(function(result) {
            console.log(result);
            $("#data").html(result);
        }).fail(function(jqXHR, textStatus, errorThrown) {
            console.log("fail: ",textStatus, errorThrown);
        });
    });
});



