import sys
sys.path.append("../")
sys.path.append("./DeLinker/analysis/")
sys.path.append("./DeLinker/")

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem.Draw import MolDrawing, DrawingOptions, MolToFile
from rdkit.Chem import MolStandardize

import numpy as np

from itertools import product
from joblib import Parallel, delayed
import re
from collections import defaultdict

from IPython.display import clear_output
IPythonConsole.ipython_useSVG = True

from DeLinker_test import DenseGGNNChemModel
from analysis import frag_utils, rdkit_conf_parallel
from data.prepare_data import read_file, preprocess
from examples import example_utils

import glob
import pandas as pd
import json
import os




def run_DeLinker(webapp_session_dir, coords_to_replace, 
                         data_path_fname = "./DeLinker_test_data",
                         json_out_fname = "DeLinker_test", 
                         num_elabs = 250,
                         output_smi_id = 'gen_linkers', image_fname = 'linked_mols.svg'):
    
    
    #Load in molecule(s)
    mols = [Chem.MolFromMolFile(x) for x in glob.glob(f'{webapp_session_dir}/ligand*.sdf')]
    
    #Load in user specified inputs
    if not os.path.exists(f'{webapp_session_dir}/user_input.json'):
        #i.e. the user hasn't specified any inputs:
        if len(mols) > 1:
            return 'Must specify a linker length if more than one fragment is provided'
        elif len(mols) == 1:
            num_elabs = 250
            linker_length = coords_to_replace.shape[0]
        
    else:
        #i.e. there is a user input file
        with open(f'{webapp_session_dir}/user_input.json', 'r') as f:
            user_input = json.load(f)
        
        if len(mols) > 1 and "max_linker_length" not in user_input.keys():
            return 'Must specify a linker length if more than one fragment is provided'
        elif len(mols) == 1 and "max_linker_length" not in user_input.keys():
            linker_length = coords_to_replace.shape[0]
        else:
            linker_length = int(user_input["max_linker_length"])
            
        if "num_linkers" not in user_input.keys():
            print("Assuming default number of elaborations")
        else:
            num_elabs = int(user_input["num_linkers"])
        
        
    
   
    
  
    #Obtain substructure we're trying to replace via the atomic coordinates
    linker, fragments, mol_smiles, mol_out = match_coords_to_indices(mols, coords_to_replace)
    

    
    #Prepare the model to be run
    model_args = process_pharmacophoric_info(mol_smiles, linker, fragments, mol_out,
                                           data_path_fname=f'{webapp_session_dir}/DeLinker_test_data.smi',
                                           json_out_fname=f'{webapp_session_dir}/molecules_DeLinker_test',
                                           output_smi_file=f'{webapp_session_dir}/DeLinker_{output_smi_id}.smi',
                                           len_linker=linker_length,
                                           num_elabs=num_elabs
                                           )
  
    #Make elabs
    generate_linkers(model_args)

    #Summary Statistics
    mol_img = summary_stats(output_smi_file=f'{webapp_session_dir}/DeLinker_{output_smi_id}.smi', image_fname=f'{webapp_session_dir}/{image_fname}')

    return mol_img


def check_DeLinker_input(webapp_session_dir, coords_to_replace, 
                             image_fname = 'DeLinker_input.svg'):
    
    
    #print(webapp_session_dir, coords_to_replace)
    #Load in molecule(s)
    mols = [Chem.MolFromMolFile(x) for x in glob.glob(f'.{webapp_session_dir}/ligand*.sdf')]
    #mol = Chem.MolFromMolFile(sdf_file)
    
    
    #Obtain substructure we're trying to replace via the atomic coordinates
    linker, fragments, mol_smiles, mol_out = match_coords_to_indices(mols, coords_to_replace)
    
    print(fragments)
    
    mol_out_2d = Chem.MolFromSmiles(fragments)
    MolToFile(mol_out_2d, f'.{webapp_session_dir}/{image_fname}', size = (500, 500))

    #Summary Statistics
    #mol_img = summary_stats(output_smi_file=f'{webapp_session_dir}/DeLinker_{output_smi_id}.smi', image_fname=f'{webapp_session_dir}/{image_fname}')

    return 0







def generate_linkers(args):
    model = DenseGGNNChemModel(args)
    model.train()
  
    model = '' #Free up some memory

    return 0


def process_pharmacophoric_info(full_smi, linker_smi, fragments_smi, original_mol,
                                data_path_fname,
                                json_out_fname, 
                                output_smi_file,
                                len_linker,
                                num_elabs = 250):
        
        
    dist, ang = frag_utils.compute_distance_and_angle(original_mol, linker_smi, fragments_smi)
        
    with open(data_path_fname, 'w') as f:
        f.write("%s %s %s" % (fragments_smi, dist, ang))
        
    raw_data = read_file(data_path_fname)
    preprocess(raw_data, "zinc", None, True, full_name = json_out_fname)    
    
    # Arguments for DeLinker
    args = defaultdict(None)
    args['--dataset'] = 'zinc'
    args['--config'] = '{"generation": true, \
                         "batch_size": 1, \
                         "number_of_generation_per_valid": %d, \
                         "min_atoms": %d, "max_atoms": %d, \
                         "train_file": "%s.json", \
                         "valid_file": "%s.json", \
                         "output_name": "%s"}' % (num_elabs, len_linker, len_linker, 
                                                              json_out_fname, json_out_fname, 
                                                              output_smi_file)
    args['--freeze-graph-model'] = False
    args['--restore'] = './DeLinker/models/pretrained_DeLinker_model.pickle'
    

    return args

def match_coords_to_indices(mols, coords, linking = True):
  
    if len(mols) == 1:
        #i.e. only a single molecule provided and 
        #the user has clicked a substructure to be removed
        
    
        
        mol = mols[0]
        
        print(Chem.MolToSmiles(mol))
        
        remove_indices = []
        for c in coords:
    
            for atom in mol.GetAtoms():
                if np.linalg.norm(np.array(mol.GetConformer().GetAtomPosition(atom.GetIdx())) - c) < 0.05:
                    remove_indices.append(atom.GetIdx())
    
        print(f'Indices to be removed from Molecule: {remove_indices}')
        elab_ev = []
        frag_ev = []
        for idx in remove_indices:
            for nei in mol.GetAtomWithIdx(idx).GetNeighbors():
                if nei.GetIdx() not in remove_indices:
                    elab_ev.append(idx)
                    frag_ev.append(nei.GetIdx())
    
        for idx in range(len(elab_ev)):
            print(f'Exit vector pair {idx + 1}: {frag_ev[idx]}, {elab_ev[idx]}')
                    
    
        bonds_to_break = [mol.GetBondBetweenAtoms(elab_ev[x],frag_ev[x]).GetIdx() for x in range(len(elab_ev))]
        fragmented_mol = Chem.FragmentOnBonds(mol, bonds_to_break)
        fragmentation = Chem.MolToSmiles(fragmented_mol).split('.')
        
        if linking:
            fragments = []
            for fragment in fragmentation:
                if len([x for x in fragment if x =="*"]) ==2:
                    linker=fragment
                else:
                    fragments.append(fragment)
            fragments = '.'.join(fragments)
            linker = re.sub('[0-9]+\*', '*', linker)
            fragments = re.sub('[0-9]+\*', '*', fragments) 
        
            return linker, fragments, Chem.MolToSmiles(mol), mol
        
    elif len(mols) == 2:
        #i.e. two fragments have been provided and we just want to link them
        
      
        
        if coords.shape[0] != 2:
            raise ValueError('For linking two separate fragments you can only select two atoms!')
         
            
        print(Chem.MolToSmiles(mols[0]))
        print(Chem.MolToSmiles(mols[1]))
        
        combo_no_exit = Chem.CombineMols(mols[0],mols[1])
        combo = Chem.CombineMols(combo_no_exit, Chem.MolFromSmiles("*.*"))
        edcombo = Chem.EditableMol(combo)
        num_heavy_atoms = combo.GetNumHeavyAtoms()
        
        add_dummy_atoms = []
        
        for atom in combo.GetAtoms():
            for c in coords:
                if np.linalg.norm(np.array(combo.GetConformer().GetAtomPosition(atom.GetIdx())) - c) < 0.05:
                    add_dummy_atoms.append(atom.GetIdx())
        
        
        edcombo.AddBond(num_heavy_atoms, add_dummy_atoms[0], order=Chem.rdchem.BondType.SINGLE)
        edcombo.AddBond(num_heavy_atoms+1, add_dummy_atoms[1], order=Chem.rdchem.BondType.SINGLE)
        editedcombo = edcombo.GetMol()
        _ = AllChem.Compute2DCoords(editedcombo)
        Chem.SanitizeMol(editedcombo)
        
        print(Chem.MolToSmiles(combo))
        
        mol_to_link = edcombo.GetMol()
        Chem.SanitizeMol(mol_to_link)
        
        # Convert exit vectors to carbons for conformer generation
        du = Chem.MolFromSmiles('*')
        mol_to_link_carbon = AllChem.ReplaceSubstructs(mol_to_link,du,Chem.MolFromSmiles('C'),True)[0]
        Chem.SanitizeMol(mol_to_link_carbon)
        # Generate conformer
        mol_to_link_carbon = Chem.AddHs(mol_to_link_carbon)
        AllChem.ConstrainedEmbed(mol_to_link_carbon, combo_no_exit, randomseed=42)
        mol_to_link_carbon = Chem.RemoveHs(mol_to_link_carbon)


        # Add this conformer to the two unlinked fragments
        conf = mol_to_link.GetConformer()
        ref_conf = mol_to_link_carbon.GetConformer()
        for i in range(mol_to_link_carbon.GetNumAtoms()):
            pos = list(ref_conf.GetAtomPosition(i))
            conf.SetAtomPosition(i, pos)
        conf.SetId(0)
        _ = mol_to_link.AddConformer(conf)



        return "", Chem.MolToSmiles(mol_to_link), "", mol_to_link
        
        
        


def summary_stats(output_smi_file, image_fname = None):

    DeLinker_mols = pd.read_csv(output_smi_file, sep = ' ')
    DeLinker_mols.columns = ['frag', 'full', 'gen']

    gen_vc = DeLinker_mols['gen'].value_counts()
    #img = Chem.Draw.MolsToGridImage([Chem.MolFromSmiles(x) for x in list(gen_vc.index)[:12]], molsPerRow = 3, returnPNG=False)
    
    #img.save(image_fname)
    mols_per_row=3

    svg = Draw.MolsToGridImage([Chem.MolFromSmiles(x) for x in list(gen_vc.index)[:12]], molsPerRow=mols_per_row, useSVG=True)
    #svg.save(image_fname)
    with open(image_fname, 'w') as f:
        f.write(svg.data)

    #if filename is not None:
    #    if not filename.endswith('.svg'):
    #        filename += '.svg'
    #    with open(filename, 'w') as f:
    #        f.write(svg)
    return svg 

    #return img
