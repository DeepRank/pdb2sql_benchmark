import numpy as np 
import subprocess as sp
from pdb2sql.StructureSimilarity import StructureSimilarity
import matplotlib.pyplot as plt
import time
import os 


case = 'BL00130001'
# case = 'BL00070001'

decoy = case + '_decoy.pdb'
ref =  case + '_ref.pdb'


deep_data = {}


print('\n-->' + decoy)

sim = StructureSimilarity(decoy,ref, enforce_residue_matching=False)
lrmsd = sim.compute_lrmsd_fast(method='svd',lzone='1KLG.lzone')
irmsd = sim.compute_irmsd_fast(method='svd',izone='1KLG.izone')
fnat = sim.compute_fnat_fast()

dockQ = sim.compute_DockQScore(fnat,lrmsd,irmsd)


mol_name = os.path.basename(decoy).split('.')[0]
deep_data[mol_name] =  [fnat,lrmsd,irmsd]



print("\n pdb2sql    : fnat = %1.6f\tlrmsd = %2.7f\tirmsd = %2.7f" %(deep_data[mol_name][0],deep_data[mol_name][1],deep_data[mol_name][2]))


	




