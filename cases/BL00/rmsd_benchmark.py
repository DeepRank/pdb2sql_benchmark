import numpy as np 
import subprocess as sp
from pdb2sql.StructureSimilarity import StructureSimilarity
# import matplotlib.pyplot as plt
import time
import os 




decoys = 'decoys'
ref = os.path.join('ref', 'BL00190001.pdb')
decoy_list = [ os.path.join(decoys,d) for d in  os.listdir(decoys) if d.endswith('.pdb')]



deep_data = {}
t0 = time.time()
for i,decoy in enumerate(decoy_list):
	print('\n-->' + decoy)


	sim = StructureSimilarity(decoy,ref, enforce_residue_matching=False)
	lrmsd_fast = sim.compute_lrmsd_fast(method='svd',lzone='BLO.lzone') #, name=['CA'])
	lrmsd = sim.compute_lrmsd_pdb2sql()#name=['CA'])
	irmsd_fast = sim.compute_irmsd_fast(method='svd',izone='BLO.izone')
	irmsd = sim.compute_irmsd_pdb2sql()
	# fnat = sim.compute_fnat_fast()

	# dockQ = sim.compute_DockQScore(fnat,lrmsd,irmsd)

	mol_name = os.path.basename(decoy).split('.')[0]
	deep_data[mol_name] =  [lrmsd_fast,lrmsd, irmsd_fast, irmsd]


	print("\nlrmsd_pdb2sql = %2.7f\tlrmsd_fast = %2.7f\nirmsd_pdb2sql = %2.7f\tirmsd_fast = %2.7f" %(deep_data[mol_name][0],deep_data[mol_name][1],deep_data[mol_name][2],deep_data[mol_name][3]))


t1=time.time()-t0
print('total time %f' %t1 )

	




