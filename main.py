#!/usr/bin/env python
# coding: utf-8

# In[1]:
import re,sys,os,getopt
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.six import StringIO
import pandas as pd
import openbabel
from rdkit.Chem import PandasTools
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Descriptors
import generate_fp_test as fpg
from mpds_id_generator import mpds_id_gen

# Read SMILES strings from a text file
smiles_list = []
with open('dataset.txt', 'r') as f:
    smiles_list = f.readlines()

# Iterate through the SMILES strings and calculate molecular weights
for smiles in smiles_list:
    try:
        mol = Chem.MolFromSmiles(smiles)  # Convert SMILES to RDKit molecule
        if mol is not None:
            mol_wt = Descriptors.MolWt(mol)  # Calculate molecular weight
            op = fpg.fp_gen(smiles,mol_wt)
            mpds_id_gen()
            print(f"Molecular Weight for SMILES: {mol_wt:.2f}")
        else:
            print(f"Error: Invalid SMILES string at index")
    except Exception as e:
        print(f"An error occurred for SMILES: {e}")
print('<-- The program completed successfully for your given dataset -->')