import re,sys,os,getopt,openbabel
from openbabel import pybel
import rdkit as rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
import checknature as ck
import ringcount as rc
import cleansmiles as cs
import check_for_frag as ckf
r_linker_r=pybel.Smarts("[R][!R][R]")
#Patterns for mono classes
def checkpo(bis,bis2,o1,smiles,fp_list,mymol,m,fusedpat,fusedpat_res,conpat,conpat_res):
	ringcount=rc.rc(smiles)
	if (len(bis)>0):
		bs=[]
		labels=[]
		for bi in bis:
			b = m.GetBondBetweenAtoms(bi[0],bi[1])
			if b.GetBeginAtomIdx()==bi[0]:
				labels.append((10,1))
			else:
				labels.append((1,10))
			bs.append(b.GetIdx())
		nm = Chem.FragmentOnBonds(m,bs,dummyLabels=labels)
		frag = Chem.MolToSmiles(nm,True)
#Fragmented frags
		fment=frag.split(".")
		nRings=[]
		for fragment in fment:
			corrected_smiles=cs.remove_non_alphabetic(fragment)
			largest=0
			fragmentcount=rc.rc(corrected_smiles)
			if fragmentcount > largest:
				if fragmentcount>=5:
					fp_list[50]='1'
	elif len(bis2)==0:
		if ringcount>=5:
			fp_list[50]='1'
	elif ringcount>=5:
		fp_list[50]='1'
	else:
		pass
	return fp_list
