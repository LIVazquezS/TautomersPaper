'''
@author: Oliver Unke
'''

#!/usr/bin/env python3
from openbabel import openbabel as ob
ob.obErrorLog.SetOutputLevel(0) #suppress non-critical warnings from openbabel (Failed to kekulize...)
from bottom_up_amon_generator import bottom_up_amons_of
import numpy as np
import pandas as pd
#input smiles
PC9 = np.loadtxt('smiles_PC92.txt', dtype='U',comments='!')

#generate converter
obConversion = ob.OBConversion()
obConversion.SetInAndOutFormats("smi", "can") #input smiles, output canonical smiles

moles_PC9 = []
#generate OBMol object from smiles
for smiles in PC9: 
    mol = ob.OBMol() 
    obConversion.ReadString(mol, smiles)
    moles_PC9.append(mol)

def get_amons(mol,n): 
    amons = [] 
    for amon in bottom_up_amons_of(mol, n): #iterator that gives back OBMol objects for amons up to 8 heavy atoms
        smiles = obConversion.WriteString(amon).strip() #convert OBMol to canonical smiles for printing
        amons.append(smiles)
    return amons

amons_PC9 = []
for i in range(0,len(moles_PC9)):
    x = get_amons(moles_PC9[i],7)
    amons_PC9.append(x)
    
dictionary = {'smiles':PC9, 'amons':amons_PC9}
df = pd.DataFrame(dictionary)
df.to_csv('Amons_PC9.csv')
