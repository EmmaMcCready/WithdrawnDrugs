#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 19 10:45:38 2021

@author: emmamcc
"""

# First, relevant packages need to be imported:
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Crippen
from rdkit.Chem import Descriptors


# Read in the dataset. The .csv file is separated by a semicolon.
file = pd.read_csv('/home/emmamcc/Documents/WithdrawnDrugs/withdrawn_dataset.csv', sep=";")


# For convenience, separate the column data into  lists
smiles = file["SMILES"].tolist()


# Assigning data frame variable for each descriptor:
mol=[]
h_acc=[]
h_don=[]
totalPSA=[]
rot=[]
logP=[]

for smiles in smiles:
        # using rdkit to read in the smiles
        chem_smiles = Chem.MolFromSmiles(smiles)
    
        # Calculating molecular weight using descriptors and appending into a list
        mol_wt = Chem.Descriptors.ExactMolWt(chem_smiles)
        mol.append(mol_wt)
    
        # Calculating Number of hydrogen acceptors using descriptors and appending into a list
        NumHAcceptors = Chem.Descriptors.NumHAcceptors(chem_smiles)
        h_acc.append(NumHAcceptors)
        
        # Calculating Number of hydrogen donors using descriptors and appending into a list
        NumHDonors = Chem.Descriptors.NumHDonors(chem_smiles)
        h_don.append(NumHDonors)
        
        # Calculating total polar surface area (tPSA) and appending into a list
        totalpolarsurfacearea = Chem.rdMolDescriptors.CalcTPSA(chem_smiles)
        totalPSA.append(totalpolarsurfacearea)
        
        # Calculating Number of rotatable bonds using descriptors and appending into a list
        NumRotatableBonds= Chem.Descriptors.NumRotatableBonds(chem_smiles)
        rot.append(NumRotatableBonds)
        
        # Calculating logP using Crippen's approach, reference: Wildman, S. and G. Crippen. “Prediction of Physicochemical Parameters by Atomic Contributions.” J. Chem. Inf. Comput. Sci. 39 (1999): 868-873.
        LOGP = Chem.Crippen.MolLogP(chem_smiles)
        logP.append(LOGP)
        


# Creating a table using pandas
physiochem_properties = pd.DataFrame(list(zip(mol,h_acc,h_don,rot,totalPSA,logP)), columns = ['Mol_wt','Num_H_Acceptors','Num_H_Donors','Num_Rotatable_Bonds','Total_PSA','LogP'])
print(physiochem_properties)


# Lastly, once I was sure I was happy with the data, I appended this dataframe to the list of drugs
result = pd.concat([file, physiochem_properties], axis=1)

result.to_csv("withdrawn_dataset_updated.csv", sep=';')