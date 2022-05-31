from flask import Flask, request, render_template, jsonify
from flask_cors import CORS
from rdkit import Chem
from rdkit.Chem.Draw import MolToFile
import shutil
import numpy as np

app = Flask(__name__)
CORS(app)

app.config['UPLOAD_FOLDER'] = './static'

@app.route('/', methods= ['GET', 'POST'])
def get_message():
    # if request.method == "GET":
    print("Got request in main function")
    return render_template("index.html")

@app.route('/upload_static_file', methods=['POST'])
def upload_static_file():
    print("Got request in static files") 
    print(request.files)
    f = request.files['static_file']
    print(f.filename)
    f_out = f'static/imgs/{f.filename}'
    f.save(f_out) 
    
    shutil.copyfile(f_out, f'static/imgs/ligand.sdf')
    
    
    mol = import_mol(f_out)
    get_mol_image_with_atom_numbers(mol, f'static/imgs/mol_with_indices.png')
    
    resp = {"success": True, "response": "file saved!"}
    return jsonify(resp), 200


@app.route('/save_array',methods=["GET", "POST"])
def save_array():
    if request.method == 'POST':
        data = request.json
        print('Printing Json...')
        print(data)

        for idx, coords in enumerate(data):
            data[idx] = [float(x) for x in coords]

        print(np.array(data))
        
    
        #Print norms:
        for coords in data:
            print(np.linalg.norm(coords))
        
        
        #with open('static/imgs/coords.txt', 'w') as f:
            #f.write(data)
            


        return jsonify(data)
    return render_template("index.html")



def import_mol(fpath):
    return Chem.MolFromMolFile(fpath)

def do_canonical_atom_renumbering(mol):
    #Get canonical atom ordering and renumber 
    #return the renumbered molecule
    mol_neworder = tuple(zip(*sorted([(j, i) for i, j in enumerate(Chem.CanonicalRankAtoms(mol))])))[1]#
    mol_renum = Chem.RenumberAtoms(mol, mol_neworder)
    
    return mol_renum

def add_atom_indices(mol):
    for i, a in enumerate(mol.GetAtoms()):
        a.SetAtomMapNum(i)


def get_mol_image_with_atom_numbers(mol, img_file_name):
    
    mol = Chem.MolFromSmiles(Chem.MolToSmiles(mol))
    mol = do_canonical_atom_renumbering(mol)
    add_atom_indices(mol)
    MolToFile(mol, img_file_name, size = (500, 500))
    return 0





if __name__ == "__main__":
    app.run(host='0.0.0.0', debug=True)
