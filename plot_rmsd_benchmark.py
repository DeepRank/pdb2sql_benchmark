import numpy as np 
import subprocess as sp
from pdb2sql.StructureSimilarity import StructureSimilarity
import matplotlib.pyplot as plt
import time




deep = np.loadtxt('deep.dat')
hdk = np.loadtxt('hdk.dat')


plt.scatter(hdk[:,2],deep[:,2],label='i-rmsd')
mini = np.min(deep[:,2])
maxi = np.max(deep[:,2])
plt.plot(  [mini,maxi],[mini,maxi],'--',color='black' )
plt.grid()
plt.legend(loc=4)
plt.xlabel('ProFit')
plt.ylabel('pdb2sql')
plt.show()




