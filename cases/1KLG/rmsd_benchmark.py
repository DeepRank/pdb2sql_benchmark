import numpy as np 
import subprocess as sp
from pdb2sql.StructureSimilarity import StructureSimilarity
import matplotlib.pyplot as plt
import time
import os 

decoys = 'decoys'
ref = os.path.join('ref', '1KLG.pdb')
decoy_list = [ os.path.join(decoys,d) for d in  os.listdir(decoys) if d.endswith('.pdb')]

deep_data = {}
t0 = time.time()
for i,decoy in enumerate(decoy_list):
	print('\n-->' + decoy)

	sim = StructureSimilarity(decoy,ref, enforce_residue_matching=False)
	lrmsd = sim.compute_lrmsd_fast(method='svd',lzone='1KLG.lzone')
	irmsd = sim.compute_irmsd_fast(method='svd',izone='1KLG.izone')
	fnat = sim.compute_fnat_fast()#ref_pairs='1AK4.refpaircd s')
	#fnat = sim.compute_Fnat_pdb2sql()
	dockQ = sim.compute_DockQScore(fnat,lrmsd,irmsd)

	# mol_name = decoy.split('/')[-1].split('.')[0]
	mol_name = os.path.basename(decoy).split('.')[0]
	deep_data[mol_name] =  [fnat,lrmsd,irmsd]
	np.savetxt(mol_name+'.LRMSD',[lrmsd])
	np.savetxt(mol_name+'.IRMSD',[irmsd])
	np.savetxt(mol_name+'.FNAT',[fnat])
	np.savetxt(mol_name+'.DOCKQ',[dockQ])

	# deep[i,:] = deep_data[mol_name]
	# hdk[i,:] = haddock_data[mol_name]

	# print("HADDOCK : fnat = %1.6f\tlrmsd = %2.7f\tirmsd = %2.7f" %(haddock_data[mol_name][0],haddock_data[mol_name][1],haddock_data[mol_name][2]))
	print("DEEP    : fnat = %1.6f\tlrmsd = %2.7f\tirmsd = %2.7f" %(deep_data[mol_name][0],deep_data[mol_name][1],deep_data[mol_name][2]))
	# print("DOCKQ   : %f" %dockQ)

t1=time.time()-t0
print('total time %f' %t1 )

	




