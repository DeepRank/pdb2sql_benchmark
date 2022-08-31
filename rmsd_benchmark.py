import numpy as np 
import subprocess as sp
from pdb2sql.StructureSimilarity import StructureSimilarity
import matplotlib.pyplot as plt
import time

'''
WHERE TO FIND HADDOCK DATA 

L-RMSD
ON ALCAZAR
/home/benchmark/docking-benchmark4/runs-cmrestraints/<MOL NAME>/run1/structures/it1/water/l-RMSD.dat

FNAT
LOCAL
~/Documents/projects/deeprank/data/HADDOCK/BM4_dimers/model_qualities/Fnat/water/<MOL>.Fnat

IRMSD
LOCAL
Documents/projects/deeprank/data/HADDOCK/BM4_dimers/model_qualities/i-rmsd/water/<MOL>.irmsd
'''

BM4 = '/home/nico/Documents/projects/deeprank/data/HADDOCK/BM4_dimers/' 
decoys = BM4 + '/decoys_pdbFLs/1AK4/water/'
ref = BM4 + '/BM4_dimers_bound/pdbFLs_ori/1AK4.pdb'


decoy_list = sp.check_output('ls %s/*.pdb' %decoys,shell=True).decode('utf-8').split()


haddock_data = {}
haddock_files =  ['1AK4.Fnat','1AK4.lrmsd','1AK4.irmsd']


for i,fname in enumerate(haddock_files):

	f = open(fname,'r')
	data = f.readlines()
	f.close()

	for line in data:

		if line[0] == '#':
			continue

		line = line.split()
		mol_name = line[0].split('.')[0]

		if i == 0:
			haddock_data[mol_name] = np.zeros(3)

		haddock_data[mol_name][i] = float(line[1])


nconf = len(haddock_data)
deep = np.zeros((nconf,3))
hdk = np.zeros((nconf,3))

deep_data = {}
t0 = time.time()
for i,decoy in enumerate(decoy_list):

	print('\n-->' + decoy)

	sim = StructureSimilarity(decoy,ref)
	lrmsd = sim.compute_lrmsd_fast(method='svd',lzone='1AK4.lzone')
	irmsd = sim.compute_irmsd_fast(method='svd',izone='1AK4.izone')
	fnat = sim.compute_fnat_fast()#ref_pairs='1AK4.refpairs')
	#fnat = sim.compute_Fnat_pdb2sql()
	dockQ = sim.compute_DockQScore(fnat,lrmsd,irmsd)

	mol_name = decoy.split('/')[-1].split('.')[0]

	deep_data[mol_name] =  [fnat,lrmsd,irmsd]
	np.savetxt(mol_name+'.LRMSD',[lrmsd])
	np.savetxt(mol_name+'.IRMSD',[irmsd])
	np.savetxt(mol_name+'.FNAT',[fnat])
	np.savetxt(mol_name+'.DOCKQ',[dockQ])

	deep[i,:] = deep_data[mol_name]
	hdk[i,:] = haddock_data[mol_name]

	print("HADDOCK : fnat = %1.6f\tlrmsd = %2.7f\tirmsd = %2.7f" %(haddock_data[mol_name][0],haddock_data[mol_name][1],haddock_data[mol_name][2]))
	print("DEEP    : fnat = %1.6f\tlrmsd = %2.7f\tirmsd = %2.7f" %(deep_data[mol_name][0],deep_data[mol_name][1],deep_data[mol_name][2]))
	print("DOCKQ   : %f" %dockQ)

t1=time.time()-t0
print('total time %f' %t1 )

	


np.savetxt('deep.dat',deep)
np.savetxt('hdk.dat',hdk)


plt.subplot(3,1,1)
plt.scatter(hdk[:,0],deep[:,0],label='Fnat')
mini = np.min(deep[:,0])
maxi = np.max(deep[:,0])
plt.plot(  [mini,maxi],[mini,maxi],'--',color='black' )
plt.legend(loc=4)
plt.xlabel('PROFIT')
plt.ylabel('DEEP')

plt.subplot(3,1,2)
plt.scatter(hdk[:,1],deep[:,1],label='l-rmsd')
mini = np.min(deep[:,1])
maxi = np.max(deep[:,1])
plt.plot(  [mini,maxi],[mini,maxi],'--',color='black' )  
plt.legend(loc=4)
plt.xlabel('PROFIT')
plt.ylabel('DEEP')

plt.subplot(3,1,3)
plt.scatter(hdk[:,2],deep[:,2],label='i-rmsd')
mini = np.min(deep[:,2])
maxi = np.max(deep[:,2])
plt.plot(  [mini,maxi],[mini,maxi],'--',color='black' )
plt.legend(loc=4)
plt.xlabel('PROFIT')
plt.ylabel('DEEP')
plt.show()




