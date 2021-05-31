#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  10 15:24:21 2020

@author: L.I.Vazquez-Salazar
@email: litzavazquezs@gmail.com

Example of a file for the calculation and plotting of KL divergency.
"""
import numpy as np
import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
from metrics import Metrics

#Read all the files related to the bonds with C.
bonds = []
for path, dirs, files in os.walk("."):
    for file in files:
        if file.startswith('bonds_C') and file.endswith(".csv"):
            bonds.append(path+"/"+file)


def read_files(array):
    '''
    Read the list of files and separate them by target or test.

    Parameters
    ----------
    array : List of files

    Returns
    -------
    values : Dictionary that refeers to each of the dataframes read it. 

    '''
    values = {}
    for item in array:
        if item.find('QM9.csv') != -1:
            values['QM9'] = pd.read_csv(item)
        elif item.find('PC9.csv') != -1:
            values['PC9'] = pd.read_csv(item)
        elif item.find('ANIE.csv') != -1:
            values['ANIE'] = pd.read_csv(item)
        elif item.find('PC9_9p.csv') != -1:
            values['PC9 opt'] = pd.read_csv(item)
        elif item.find('QM9_9p.csv') != -1:
            values['QM9 opt'] = pd.read_csv(item)
        elif item.find('ANIE_9p.csv') != -1:
            values['ANIE opt'] = pd.read_csv(item)
            
    return values

bonds_read = read_files(bonds)
 
  
def getdata(dictionary,key):
    '''
    Separate the data for a type of bond for the different databases evaluated.
    This considers only the reference distribution.
    

    Parameters
    ----------
    dictionary : Created by read_files. It contains the df for each of the databases.
    key : It is name of the bond on the databases

    Returns
    -------
    df : new database for the type of bond selected.

    '''
    X_ANIE = dictionary['ANIE'][key]
    X_PC9 = dictionary['PC9'][key]
    X_QM9 = dictionary['QM9'][key]

    dictionary = {'QM9':X_QM9, 'PC9':X_PC9, 'ANIE':X_ANIE}
    df = pd.DataFrame(dictionary)
    return df

def getdata_test(dictionary,key):
    '''
    Separate the data for a type of bond for the different databases evaluated.
    This considers only the target distribution.
    

    Parameters
    ----------
    dictionary : Created by read_files. It contains the df for each of the databases.
    key : It is name of the bond on the databases

    Returns
    -------
    df : new database for the type of bond selected.

    '''   
    X_Tauto_opt_ANIE = dictionary['ANIE opt'][key]
    X_Tauto_opt_PC9 = dictionary['PC9 opt'][key]
    X_Tauto_opt_QM9 = dictionary['QM9 opt'][key]
    
    dictionary = {'Test set opt QM9':X_Tauto_opt_QM9, 'Test set opt PC9':X_Tauto_opt_PC9
                  , 'Test set opt ANIE':X_Tauto_opt_ANIE }
    df = pd.DataFrame(dictionary)
    return df


bonds_CC = getdata(bonds_read,'CC')
test_CC = getdata_test(bonds_read,'CC')

bonds_CO = getdata(bonds_read, 'CO')
test_CO = getdata_test(bonds_read, 'CO')

bonds_CN = getdata(bonds_read, 'CN')
test_CN = getdata_test(bonds_read, 'CN')

#Values for the different types of bonds considered.
Min_CC = [1.167,1.280,1.400,1.458,1.374,1.350,1.470,1.479,1.436,1.441,1.425,1.430]
Max_CC = [1.197,1.405,1.568,1.610,1.474,1.440,1.538,1.539,1.481,1.512,1.441,1.448]

Min_CN = [1.482,1.446,1.461,1.314,1.369,1.461,1.340,1.422,1.311,1.273,1.325,1.300,1.140,1.131]
Max_CN = [1.510,1.572,1.506,1.419,1.384,1.470,1.476,1.442,1.324,1.339,1.369,1.348,1.148,1.449]

Min_CO = [1.395,1.405,1.417,1.435,1.430,1.324,1.341,1.279,1.328,1.379,1.332,1.353,1.363,1.375,
          1.394,1.188,1.232,1.203,1.181,1.184,1.187,1.193]
Max_CO = [1.449,1.458,1.438,1.501,1.501,1.342,1.363,1.320,1.420,1.393,1.377,1.373,1.377,1.391
          ,1.408,1.238,1.262,1.241,1.207,1.193,1.187,1.243]

def compute_KL(reference,target,array_min,array_max):
    '''
    Read the created databases for a specific type of bond and then compute 
    the value of the KL divergency.

    Parameters
    ----------
    reference : Dataframe of the values of the bond lenghts of a specific type of bond
                on the reference databases.
    target : Datafram of the values of the bond lenghts of a specific type of bond
                on the target database.

    Returns
    -------
    X_values : Evaluated values.
    KL_values : Values of the cummulative KL

    '''
    mt = Metrics(reference,target)
    data_QM9 = mt.get_data('QM9','Test set opt QM9')
    data_PC9 = mt.get_data('PC9','Test set opt PC9')
    data_ANIE = mt.get_data('ANIE','Test set opt ANIE')
    
    QM9 = []
    PC9 = []
    ANIE =[]
    for i in range(0,len(array_min)):
        QM9_kl = mt.KL_divergence_int(data_QM9,array_min[i],array_max[i])
        PC9_kl = mt.KL_divergence_int(data_PC9,array_min[i],array_max[i])
        ANIE_kl = mt.KL_divergence_int(data_ANIE,array_min[i],array_max[i])
        QM9.append(QM9_kl)
        PC9.append(PC9_kl)
        ANIE.append(ANIE_kl)
    KL_values = [QM9,PC9,ANIE]
    return KL_values

cc_bonds_KL = compute_KL(bonds_CC,test_CC,Min_CC,Max_CC) 
cn_bonds_KL = compute_KL(bonds_CN,test_CN, Min_CN,Max_CN)
co_bonds_KL = compute_KL(bonds_CO,test_CO,Min_CO,Max_CO)

def create_csv(array,csv_name):
    '''
    Saves the data obtained on the previous step as a .csv file.

    Parameters
    ----------
    array : Values of the KL divergency for a type of bond.
    csv_name : Name of the output .csv file

    Returns
    -------
    None.

    '''
    dct = {'QM9':array[0], 'PC9':array[1],'ANIE':array[2]}
    df = pd.DataFrame(dct)
    df.to_csv(csv_name)
    
create_csv(cc_bonds_KL,'CC_KL_9p.csv')
create_csv(cn_bonds_KL,'CN_KL_9p.csv')
create_csv(co_bonds_KL,'CO_KL_9p.csv')


