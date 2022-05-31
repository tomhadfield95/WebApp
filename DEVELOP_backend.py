def run_DEVELOP(sdf_file, coords_to_replace, 
                          frag_sdf_fname = "scaffold_elaboration_core.sdf", 
                          pharm_sdf_fname = "scaffold_elaboration_pharmacophores.sdf", 
                          data_path_fname = "./scaffold_elaboration_test_data",
                          types_fname = "scaffold_elaboration.types",
                          json_out_fname = "scaffold_elaboration_test", 
                          num_elabs = 250,
                          output_smi_file = 'scaffold_elaboration_gen'):

  #Load in molecule
  mol = Chem.MolFromMolFile(sdf_file)
  
  #Obtain substructure we're trying to replace via the atomic coordinates
  scaffold,  elaboration, full_smiles = match_coords_to_indices(mol, coords_to_replace)

  #Prepare the model to be run
  model_args = process_pharmacophoric_info(full_smiles, elaboration, scaffold, sdf_file,
                                           frag_sdf_fname=frag_sdf_fname,
                                           pharm_sdf_fname=pharm_sdf_fname,
                                           types_fname=types_fname,
                                           json_out_fname=json_out_fname,
                                           num_elabs=num_elabs,
                                           output_smi_file=output_smi_file)
  
  #Make elabs
  generate_elaborations(model_args)

  #Summary Statistics
  mol_img = summary_stats(output_smi_file=output_smi_file)

  return mol_img


def generate_elaborations(args):
  model = DenseGGNNChemModel(args)
  model.train()
  
  model = '' #Free up some memory

  return 0

def process_pharmacophoric_info(full_smi, elab_smi, scaffold_smi, full_sdf_fname,
                                frag_sdf_fname = "scaffold_elaboration_core.sdf", 
                                pharm_sdf_fname = "scaffold_elaboration_pharmacophores.sdf", 
                                data_path_fname = "./scaffold_elaboration_test_data",
                                types_fname = "scaffold_elaboration.types",
                                json_out_fname = "scaffold_elaboration_test", 
                                num_elabs = 250,
                                output_smi_file = 'scaffold_elaboration_gen'):

  # Write data to file
  with open(data_path_fname, 'w') as f:
      f.write("%s %s %s" % (full_smi, elab_smi, scaffold_smi))

  raw_data = read_file(data_path_fname, add_idx=True, calc_pharm_counts=True)
  preprocess(raw_data, "zinc", json_out_fname, "./", False)


  fragmentations_pharm, fails = frag_utils.create_frags_pharma_sdf_dataset([[full_smi, elab_smi, scaffold_smi, 0, 0]], 
                                                                          full_sdf_fname, dataset="CASF",
                                                                          sdffile=core_path,
                                                                          sdffile_pharm=pharmacophores_path,
                                                                          prot="", verbose=True)


  # Write .types file
  with open(types_fname, 'w') as f:
    f.write('1 ' + frag_sdf_fname + ' ' + pharm_sdf_fname)


  #Set up DEVELOP arguments
  args = defaultdict(None)

  args['--dataset'] = 'zinc'
  args['--config'] = '{"generation": true, \
                      "batch_size": 1, \
                      "number_of_generation_per_valid": %d, \
                      "train_file": "./molecules_%s.json", \
                      "valid_file": "./molecules_%s.json", \
                      "train_struct_file": "./%s", \
                      "valid_struct_file": "./%s", \
                      "struct_data_root": "./", \
                      "output_name": "DEVELOP_%s.smi"}' % (num_elabs, json_out_fname, 
                                                           json_out_fname, types_fname,
                                                           types_fname, output_smi_file)


  args['--freeze-graph-model'] = False
  args['--restore'] = '../models/scaffold_elaboration/pretrained_DEVELOP_model.pickle'


  return args


def match_coords_to_indices(mol, coords):
  
  elab_indices = []
  for c in coords:

    for atom in mol.GetAtoms():
      if np.linalg.norm(np.array(mol.GetConformer().GetAtomPosition(atom.GetIdx())) - c) < 0.05:
        elab_indices.append(atom.GetIdx())

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

  elaboration = re.sub('[0-9]+\*', '*', elaboration)
  scaffold = re.sub('[0-9]+\*', '*', scaffold)

  return scaffold, elaboration, Chem.MolToSmiles(mol)


def summary_stats(output_smi_file, image_fname = None):

  DEVELOP_mols = pd.read_csv(f"DEVELOP_{output_smi_file}.smi", sep = ' ')
  DEVELOP_mols.columns = ['frag', 'full', 'gen']

  gen_vc = DEVELOP_mols['gen'].value_counts()
  img = Chem.Draw.MolsToGridImage([Chem.MolFromSmiles(x) for x in list(gen_vc.index)[:12]], molsPerRow = 3)

  return img
