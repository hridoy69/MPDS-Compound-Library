**Molecular Property Diagnostic Suite – Compound Library (MPDS-CL)**


Molecular Property Diagnostic Suite-Compound Library (MPDS-CL), is an open-source galaxy-based cheminformatics web-portal which presents a structure-based classification of the molecules. A structure-based classification of nearly 150 million unique compounds has been done and deposited in MPDS-CL. A 56-bit fingerprint-based classification algorithm which led to a formation of 56 structurally well-defined classes. This is the official open source repository for the following paper:


Lijo John, Selvaraman Nagamani, Hridoy Jyoti Mahanta, S. Vaikundamani, Nandan Kumar, Asheesh Kumar, Esther Jamir, Lipsa Priyadarsinee, G. Narahari Sastry (2023): Molecular Property Diagnostic Suite Compound Library (MPDS-CL): A Structure based Classification of the Chemical Space. https://www.researchsquare.com/article/rs-3236523/latest.pdf 

**Requirements**

•	Python 3.9 & above

•	RDkit - 2022.03

•	Openbabel - v.3.0.0

•	Install these two packages in conda environment and activate which running the script. 


**Usage**

Required scripts are given in the scripts folder. 


Users can include the list of SMILES in “dataset.txt” and run “main.py” script. The output will be saved in “fp.txt” and “final_mpds_class_ouput.txt”. The “fp.txt” contains 56 bit fingerprint for each input compound. Each bit represents the presence of 56 features or class (0 absent; 1 present) in the compound. A compound can have multiple features but will belong to only 1 class based on priority. The output contains four columns

1st column – Identified MPDS class number

2nd column – Generated 56 bit FP

3rd column – Molecular weight of the molecule

4th column  - input SMILES


**People**

•	Lijo John

•	Selvaraman Nagamani 

•	Hridoy Jyoti Mahanta 

•	S. Vaikundamani

•	Nandan Kumar 

•	Asheesh Kumar

•	Esther Jamir

•	Lipsa Priyadarsinee

•	G. Narahari Sastry


**If you find this useful in your research, please cite**

L. John.; S. Nagamani.; H.J. Mahanta.; S. Vaikundamani.; N. Kumar.; A. Kumar.; E. Jamir.; L. Priyadarsinee.; G. N. Sastry (2023): Molecular Property Diagnostic Suite Compound Library (MPDS-CL): A Structure based Classification of the Chemical Space. DOI: 10.21203/rs.3.rs-3236523/v1 
