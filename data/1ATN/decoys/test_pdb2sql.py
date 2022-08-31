from deeprank.tools import pdb2sql 

sql = pdb2sql('1ATN_114w.pdb')
ind = sql.get_contact_atoms()
ind = ind[0]+ind[1]
xyz = sql.get('x,y,z',rowID=ind)