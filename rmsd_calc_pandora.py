#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 7  16:58:29 2022

@author: Farzaneh Parizi

USAGE: python rmsd_calc_pandora.py  <indir>  <outdir>

python rmsd_calc_pandora.py ./cases/1KLG 1KLG_PANDORA_objects 

"""

import os
import sys
import pickle
sys.path.append('/home/fmeiman/PANDORA')
import PANDORA
import csv

#def calc_lrmsd(decoy_path, ref_path, atoms):
#    sim =  StructureSimilarity(decoy_path, ref_path)
#    lrmsd = sim.compute_lrmsd_pdb2sql(exportpath=None, method='svd', name = atoms)
#
#    return lrmsd

indir=sys.argv[1]
fol = sys.argv[2]
files = []
fol_path= os.path.join(indir,fol)
for i in os.listdir(fol_path):
     if os.path.isfile(os.path.join(fol_path,i)) and 'results_' in i:
         files.append(i)

results_pkl_file= files[0]

with open(os.path.join(fol_path,results_pkl_file), 'rb') as inpkl:
       case_models = pickle.load(inpkl)

target_id = fol[:4]
case_dict = {x.model_path.split('/')[-1] : {'molpdf' : float(x.molpdf), 'object' : x} for x in case_models} #'model_object' : x,
final_scores = {}


print('Calculating RMSDs')

for model in case_dict:
    case_dir = fol_path
    model_id = model.split('.')[1]
    try:
        case_dict[model]['object'].calc_LRMSD('%s/%s.pdb' %(case_dir, target_id), atoms=['CA'])
        case_dict[model]['CA_lRMSD'] = case_dict[model]['object'].lrmsd
    except:
        raise Exception(model, case_dir)

    case_dict[model]['object'].calc_LRMSD('%s/%s.pdb' %(case_dir, target_id), atoms=['CA', 'C', 'O', 'N'])
    case_dict[model]['BB_lRMSD'] = case_dict[model]['object'].lrmsd

    case_dict[model]['object'].calc_LRMSD('%s/%s.pdb' %(case_dir, target_id), atoms=['CA', 'C', 'O', 'N', 'CB'])
    case_dict[model]['BB_CB_lRMSD'] = case_dict[model]['object'].lrmsd

#    case_dict[model]['object'].calc_LRMSD('%s/%s.pdb' %(case_dir, target_id), atoms=)
#    case_dict[model]['FA_lRMSD'] = case_dict[model]['object'].lrmsd

    del case_dict[model]['object']

header = ['Model', 'molpdf', 'CA_l-RMSD', 'BB_l-RMSD', 'BB_CB_lRMSD']

with open(os.path.join(fol_path,'rmsds_and_final_scores.tsv'), 'wt') as outfile:
    tw = csv.writer(outfile, delimiter='\t')
    tw.writerow(header)
    for key in case_dict:
        try:
            tw.writerow([key, case_dict[key]['molpdf'], case_dict[key]['CA_lRMSD'],
                          case_dict[key]['BB_lRMSD'], case_dict[key]['BB_CB_lRMSD']])
        except:
            tw.writerow([key, case_dict[key]['molpdf'], case_dict[key]['CA_lRMSD'],
                          'N/A', 'N/A', 'N/A'])

#molsort = sorted(final_scores.items(), key=lambda x:x[1][0])
with open(os.path.join(fol_path,'rmsds_topmolpdfs_'+ target_id + '.pkl'), 'wb') as dict_outpkl:
    pickle.dump(case_dict, dict_outpkl)

