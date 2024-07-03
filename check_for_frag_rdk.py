import re
def checknatureforfragrdk(smiles,rdk_trpatt3,rdk_trpatt4,rdk_trpatt5,rdk_trpatt6,rdk_trpatt7,rdk_trpatt8,rdk_fivear,rdk_sixar,rdk_conpat,rdk_fusedpat):
	pure_arom=0
	pure_alip=0
	mixed_arom_one=0
	mixed_arom_two=0
#Pure-arom
	if len(rdk_trpatt3)>=1 and len(rdk_trpatt8)==0:
		pure_arom=1
	if len(rdk_trpatt4)>=1 and len(rdk_trpatt8)==0:
		pure_arom=1
	if re.search(r'(c)',frag) and re.search('^((?!C|N|O|P|S).)*$',frag) and len(rdk_trpatt8)==0:
		pure_arom=1
	if re.search(r'(c|n|o|p|s)',frag) and re.search('^((?!C|N|O|P|S).)*$',frag):
		pure_arom=1
	if len(rdk_trpatt5)>=1 and len(rdk_trpatt8)==0:
		pure_arom=1
#Mixed-arom
	if len(rdk_trpatt3)>=1 and len(rdk_trpatt8)>=1:
#	if len(rdk_trpatt3)>=1 and len(rdk_trpatt8)>=1 and re.search('^((C|N|O|P|S).)*$',frag):
		print("found mixed in check_for_frag function")
		if len(rdk_fivear)==5 and len(rdk_sixar)==0:
			mixed_arom_one=1
		if len(rdk_fivear)==5 and len(rdk_sixar)==6:
			mixed_arom_two=1
		if len(rdk_sixar)>=12:
			mixed_arom_two=1
	if len(rdk_trpatt4)>=1 and len(rdk_trpatt8)>=1 and re.search('^((C|N|O|P|S).)*$',frag):
		print("found mixed in check_for_frag function")
		if len(rdk_fivear)==5 and len(rdk_sixar)==0:
			mixed_arom_one=1
		if len(rdk_fivear)==5 and len(rdk_sixar)==6:
			mixed_arom_two=1
#Pure-ali
	if len(rdk_trpatt3)==0 and len(rdk_trpatt4)==0 and len(result6)>=0 and len(result7)>=0 and len(rdk_trpatt8)>=1 and re.search('^((?!c|n|o|s|p).)*$',frag):
		pure_alip=1
#Pure-arom for the pure fused
	if len(rdk_trpatt3)>=1 and len(rdk_trpatt8)==0 and len(fusedpat_res)==4 and len(conpat_res)==0:
		pure_arom=1
	if len(rdk_trpatt4)>=1 and len(rdk_trpatt8)==0 and len(fusedpat_res)==4 and len(conpat_res)==0:
		pure_arom=1
	if re.search(r'(c)',frag) and re.search('^((?!C|N|O|P|S).)*$',frag) and len(rdk_trpatt8)==0 and len(fusedpat_res)==4 and len(conpat_res)==0:
		pure_arom=1
	if re.search(r'(c|n|o|p|s)',frag) and re.search('^((?!C|N|O|P|S).)*$',frag) and len(fusedpat_res)==4 and len(conpat_res)==0:
		pure_arom=1
	if len(rdk_trpatt5)>=1 and len(rdk_trpatt8)==0:
		pure_arom=1
#Mixed-arom
	if len(rdk_trpatt3)>=1 and len(rdk_trpatt8)>=1 and re.search('^((C|N|O|P|S).)*$',frag) and len(fusedpat_res)==4 and len(conpat_res)==0:
		if len(rdk_fivear)==5 and len(rdk_sixar)==4:
			mixed_arom_one=1
		if len(rdk_fivear)==2 and len(rdk_sixar)==6:
			mixed_arom_one=1
		if len(rdk_fivear)==4 and len(rdk_sixar)==12:
			mixed_arom_two=1
		if len(rdk_fivear)==5 and len(rdk_sixar)==8:
			mixed_arom_two=1
	if len(rdk_trpatt4)>=1 and len(rdk_trpatt8)>=1 and re.search('^((C|N|O|P|S).)*$',frag) and len(fusedpat_res)==4 and len(conpat_res)==0:	
		if len(rdk_fivear)==5 and len(rdk_sixar)==4:
			mixed_arom_one=1
		if len(rdk_fivear)==2 and len(rdk_sixar)==6:
			mixed_arom_one=1
		if len(rdk_fivear)==4 and len(rdk_sixar)==12:
			mixed_arom_two=1
		if len(rdk_fivear)==5 and len(rdk_sixar)==8:
			mixed_arom_two=1
#Pure-ali
	if len(rdk_trpatt3)==0 and len(rdk_trpatt4)==0 and len(result6)>=0 and len(result7)>=0 and len(rdk_trpatt8)>=1 and re.search('^((?!c|n|o|s|p).)*$',frag) and len(fusedpat_res)==4 and len(conpat_res)==0:
		pure_alip=1
	return pure_alip,pure_arom,mixed_arom_one,mixed_arom_two
