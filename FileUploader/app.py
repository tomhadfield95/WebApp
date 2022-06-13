from flask import Flask, request, render_template, jsonify
from flask_cors import CORS
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import MolToFile
import shutil
import numpy as np
import re
from generate import run_DeLinker


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
        
        
        np.savetxt('static/imgs/coords.txt', np.array(data))

        mol = Chem.MolFromMolFile('static/imgs/ligand.sdf')
        scaffold, elaboration, mol_smiles = match_coords_to_indices(mol, np.array(data), remove_dummy = True)

        print(elaboration)
        

        matching = mol.GetSubstructMatch(Chem.MolFromSmiles(elaboration))
        print(matching)
        print(Chem.MolToSmiles(mol))
        mol.RemoveAllConformers()
        
        MolToFile(mol, 'static/imgs/mol_substruct.png', size = (500, 500), highlightAtoms = matching)

        return jsonify(data)
    return render_template("index.html")





@app.route('/generate_linkers',methods=["GET", "POST"])
def generate_linkers():
    if request.method == 'POST':
        
        data = request.json

        
        coords = np.loadtxt('static/imgs/coords.txt')
        #mol = Chem.MolFromMolFile('static/imgs/ligand.sdf')
        
        run_DeLinker('static/imgs/ligand.sdf', coords)
        
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


def match_coords_to_indices(mol, coords, remove_dummy = False):
  
    elab_indices = []
    
    
    print(coords)
    print(Chem.MolToMolBlock(mol))
    
    
    for c in coords:

        for atom in mol.GetAtoms():
            if np.linalg.norm(np.array(mol.GetConformer().GetAtomPosition(atom.GetIdx())) - c) < 0.05:
                elab_indices.append(atom.GetIdx())
    
    print(elab_indices)

    elab_ev = -1
    frag_ev = -1
    for idx in elab_indices:
        for nei in mol.GetAtomWithIdx(idx).GetNeighbors():
            if nei.GetIdx() not in elab_indices:
                elab_ev = idx
                frag_ev = nei.GetIdx()

    print(f'Fragment Exit Vector Index: {frag_ev}')
    print(f'Elaboration Exit Vector Index: {elab_ev}')

    bond_to_break = mol.GetBondBetweenAtoms(elab_ev, frag_ev).GetIdx()
    fragmented_mol = Chem.FragmentOnBonds(mol, [bond_to_break])

    fragmentation = Chem.MolToSmiles(fragmented_mol).split('.')
    len_elab = coords.shape[0]
  

    if Chem.MolFromSmiles(fragmentation[0]).GetNumHeavyAtoms() == len_elab:
        elaboration = fragmentation[0]
        scaffold = fragmentation[1]
    else:
        elaboration = fragmentation[1]
        scaffold = fragmentation[0]

    if not remove_dummy:
        elaboration = re.sub('[0-9]+\*', '*', elaboration)
        scaffold = re.sub('[0-9]+\*', '*', scaffold)
    else:
        elaboration = remove_dummy_atom(elaboration)
        scaffold = remove_dummy_atom(scaffold)

    return scaffold, elaboration, Chem.MolToSmiles(mol)


def remove_dummy_atom(smiles):
    mol = Chem.MolFromSmiles(smiles)
    mol2 = AllChem.ReplaceSubstructs(mol, Chem.MolFromSmiles('*'), Chem.MolFromSmiles('[H]'), True)[0]
    mol3 = Chem.RemoveHs(mol2)

    return Chem.MolToSmiles(mol3)

if __name__ == "__main__":
    app.run(host='0.0.0.0', debug=True)
