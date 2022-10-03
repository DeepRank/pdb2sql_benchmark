import numpy as np 
import subprocess as sp
from pdb2sql.StructureSimilarity import StructureSimilarity
# import matplotlib.pyplo   t as plt
import time
import os 


case = 'BL00070001'
case = 'BL00050001'

decoy = case + '_decoy.pdb'
ref =  case + '_ref.pdb'


deep_data = {}
t0 = time.time()

print('\n-->' + decoy)

sim = StructureSimilarity(decoy,ref, enforce_residue_matching=False)
lrmsd_fast = sim.compute_lrmsd_fast(method='svd',lzone='1KLG.lzone', name=['CA'])
lrmsd = sim.compute_lrmsd_pdb2sql(name=['CA'])
irmsd_fast = sim.compute_irmsd_fast(method='svd',izone='1KLG.izone')
irmsd = sim.compute_irmsd_pdb2sql()
# fnat = sim.compute_fnat_fast()

# dockQ = sim.compute_DockQScore(fnat,lrmsd,irmsd)

mol_name = os.path.basename(decoy).split('.')[0]
deep_data[mol_name] =  [lrmsd_fast,lrmsd, irmsd_fast, irmsd]


print("\nlrmsd_pdb2sql = %2.7f\tlrmsd_fast = %2.7f\nirmsd_pdb2sql = %2.7f\tirmsd_fast = %2.7f" %(deep_data[mol_name][0],deep_data[mol_name][1],deep_data[mol_name][2],deep_data[mol_name][3]))



	




