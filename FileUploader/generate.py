import sys
sys.path.append("../")
sys.path.append("./DeLinker/analysis/")
sys.path.append("./DeLinker/")

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem.Draw import MolDrawing, DrawingOptions
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




def run_DeLinker(sdf_file, coords_to_replace, 
                         data_path_fname = "./DeLinker_test_data",
                         json_out_fname = "DeLinker_test", 
                         num_elabs = 250,
                         output_smi_file = 'gen_linkers', image_fname = './static/imgs/linked_mols.svg'):
    #Load in molecule
    mol = Chem.MolFromMolFile(sdf_file)
  
    #Obtain substructure we're trying to replace via the atomic coordinates
    linker, fragments, mol_smiles = match_coords_to_indices(mol, coords_to_replace)
    

    
    #Prepare the model to be run
    model_args = process_pharmacophoric_info(mol_smiles, linker, fragments, mol,
                                           data_path_fname=data_path_fname,
                                           json_out_fname=json_out_fname,
                                           num_elabs=num_elabs,
                                           output_smi_file=output_smi_file)
  
    #Make elabs
    generate_linkers(model_args)

    #Summary Statistics
    mol_img = summary_stats(output_smi_file=output_smi_file, image_fname=image_fname)

    return mol_img


def generate_linkers(args):
    model = DenseGGNNChemModel(args)
    model.train()
  
    model = '' #Free up some memory

    return 0


def process_pharmacophoric_info(full_smi, linker_smi, fragments_smi, original_mol,
                                data_path_fname = "./DeLinker_test_data",
                                json_out_fname = "DeLinker_test", 
                                num_elabs = 250,
                                output_smi_file = 'gen_linkers'):

    
    len_linker = Chem.MolFromSmiles(linker_smi).GetNumHeavyAtoms()
    dist, ang = frag_utils.compute_distance_and_angle(original_mol, linker_smi, fragments_smi)
    
        
    with open(data_path_fname, 'w') as f:
        f.write("%s %s %s" % (fragments_smi, dist, ang))
        
    raw_data = read_file(data_path_fname)
    preprocess(raw_data, "zinc", json_out_fname, True)    
    
    # Arguments for DeLinker
    args = defaultdict(None)
    args['--dataset'] = 'zinc'
    args['--config'] = '{"generation": true, \
                         "batch_size": 1, \
                         "number_of_generation_per_valid": %d, \
                         "min_atoms": %d, "max_atoms": %d, \
                         "train_file": "molecules_%s.json", \
                         "valid_file": "molecules_%s.json", \
                         "output_name": "DeLinker_%s.smi"}' % (num_elabs, len_linker, len_linker, 
                                                              json_out_fname, json_out_fname, 
                                                              output_smi_file)
    args['--freeze-graph-model'] = False
    args['--restore'] = './DeLinker/models/pretrained_DeLinker_model.pickle'
    

    return args

def match_coords_to_indices(mol, coords, linking = True):
  
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
    
        return linker, fragments, Chem.MolToSmiles(mol)


def summary_stats(output_smi_file, image_fname = None):

    DeLinker_mols = pd.read_csv(f"DeLinker_{output_smi_file}.smi", sep = ' ')
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
