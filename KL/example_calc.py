#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Example file of the calculations of the bond lenghts on QM9. 

Same file can be used for another databases with a .xyz format. 

Created on Mon Nov 23 16:16:33 2020

@author: L.I.Vazquez-Salazar
@email: litzavazquezs@gmail.com
"""

import numpy as np
import pandas as pd
from Mol_measure import Mol_measure 
import os
from tqdm import tqdm
import seaborn as sns


mm = Mol_measure()
# =============================================================================
# Read files

file_list = []
for path, dirs, files in os.walk("QM9"):
    for file in files:
        if file.endswith(".xyz"):
            file_list.append(path+"/"+file)
# =============================================================================

# =============================================================================
# Bonds

def calculator_bond(a,b):
    '''
    Compute the bonds between two atoms on a loop for all the molecules on the 
    database

    Parameters
    ----------
    a : Label of atom a
    b : Label of atom b

    Returns
    -------
    bond_lengths : Array of all the bond lenghts between those two atoms. 

    '''
    bond_lengths = []
    for file in tqdm(file_list):
        mol = mm.read_molecule(file)
        bond = mm.calculate_bond(a,b)
        bond_lengths.extend(bond)
    return bond_lengths


C_C_bonds = np.array(calculator_bond('C','C')) #Compute bonds between carbon atoms
C_O_bonds = np.array(calculator_bond('C','O')) #Compute bonds between carbon and oxygen atoms
C_N_bonds = np.array(calculator_bond('C','N')) #Compute bonds between carbon and nitrogen atoms.


dic_Carbon_bonds = {'CC':C_C_bonds, 'CO':C_O_bonds, 'CN':C_N_bonds}
df_c = pd.DataFrame.from_dict(dic_Carbon_bonds, orient='index').transpose()
df_c.to_csv('bonds_C.csv') #Save it as a .csv file that can be used later

N_N_bonds = np.array(calculator_bond('N','N'))#Compute bonds between nitrogen atoms
N_O_bonds = np.array(calculator_bond('N','O'))#Compute bonds between nitrogen and oxygen atoms

dic_Nitrogen_bonds = {'NN':N_N_bonds, 'NO':N_O_bonds}
df_n = pd.DataFrame.from_dict(dic_Nitrogen_bonds, orient='index').transpose()
df_n.to_csv('bonds_N.csv')#Save it as a .csv file that can be used later

