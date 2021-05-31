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

    dictionary = {'QM9':X_QM9, 'PC9':X_PC9, 'ANI-1E':X_ANIE}
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


def compute_KL(reference,target):
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
    data_ANIE = mt.get_data('ANI-1E','Test set opt ANIE')
    
    QM9_x, QM9_kl = mt.KL_divergence_cum(data_QM9)
    PC9_x, PC9_kl = mt.KL_divergence_cum(data_PC9)
    ANIE_x, ANIE_kl = mt.KL_divergence_cum(data_ANIE)
    
    KL_values = [QM9_kl,PC9_kl,ANIE_kl]
    X_values = [QM9_x[1:],PC9_x[1:], ANIE_x[1:]]
    return X_values, KL_values

X_CC, KL_CC = compute_KL(bonds_CC,test_CC)
X_CN, KL_CN = compute_KL(bonds_CN,test_CN)
X_CO, KL_CO = compute_KL(bonds_CO,test_CO)

x_value = np.linspace(1,2,num=999)
# =============================================================================
#%%
colors = ['firebrick','darkorange','deepskyblue']
colors_rgba = []
for i in colors:
    x = matplotlib.colors.to_rgba(i)
    colors_rgba.append(x)
names = ['QM9 KL', 'PC9 KL', 'ANIE KL']

fig, axs = plt.subplots(nrows=3, ncols=3, figsize=(15,15),sharex='col')
plt.rc('font', family='sans-serif', size=16)
ax1 = sns.kdeplot(data=bonds_CC, ax=axs[0][0], palette=colors_rgba, legend=True, common_norm=False)
ax1.set_title('C-C Bond distribution')
ax1.text(0.85,5,('Reference'),fontsize=35,rotation='vertical')
ax1.set_ylim(0,18)
ax2 = sns.kdeplot(data=test_CC,ax=axs[1][0],fill=False,common_norm=False,legend=False,palette=colors_rgba)
ax2.text(0.85,4.5,('Target'),fontsize=35,rotation='vertical')
ax2.set_ylim(0,18)
ax3 = axs[2][0]
for i in range(0,3):
    ax3.plot(X_CC[i],KL_CC[i], label=names[i],color=colors[i])

ax3.set_xlim(1.1,1.7)
ax3.set_ylabel(r'$D_{KL}(P||Q)$')
ax3.set_xlabel(r'r($\AA$)')
ax3.text(0.85,0.1,('KL-div'),fontsize=35,rotation='vertical')

ax4 = sns.kdeplot(data=bonds_CN, ax=axs[0][1], palette=colors_rgba, legend=True, common_norm=False,)
ax4.set_title('C-N Bond distribution')
ax4.set_ylabel('')
ax4.set_ylim(0,18)
ax5 = sns.kdeplot(data=test_CN,ax=axs[1][1],fill=False,common_norm=False,legend=False,palette=colors_rgba)
ax5.set_ylim(0,18)
ax5.set_ylabel('')
ax6 = axs[2][1]
for i in range(0,3):
    ax6.plot(X_CN[i],KL_CN[i], label=names[i],color=colors[i])
    
ax6.set_xlim(1.1,1.7)

ax6.set_xlabel(r'r($\AA$)')

ax7 = sns.kdeplot(data=bonds_CO, ax=axs[0][2], palette=colors_rgba, legend=True, common_norm=False)
ax7.set_title('C-O Bond distribution')
ax7.set_ylabel('')
ax7.set_ylim(0,18)

ax8 = sns.kdeplot(data=test_CO,ax=axs[1][2],fill=False,common_norm=False,legend=False,palette=colors_rgba)
ax8.set_ylabel('')
ax8.set_ylim(0,18)

ax9 = axs[2][2]
for i in range(0,3):
    ax9.plot(X_CO[i],KL_CO[i], label=names[i],color=colors[i])
    
ax9.set_xlim(1.1,1.7)

ax9.set_xlabel(r'r($\AA$)')

plt.subplots_adjust(hspace=0)
#plt.savefig('C-bonds-dist-KL_9p.pdf',bbox_inches='tight')
# plt.show()
