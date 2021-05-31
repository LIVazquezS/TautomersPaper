#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 15:25:07 2021

@author: L.I.Vazquez-Salazar
@email: litzavazquezs@gmail.com

Script to plot the error by number of data. It needs an extra file with the number of atoms
by element for each of the molecules. 
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.ticker import FormatStrFormatter

df_a_e = pd.read_csv('Num_of_atoms_by_element.csv')
ANIE = pd.read_csv('ANIE.csv')
PC9 = pd.read_csv('PC9.csv')
QM9 = pd.read_csv('QM9.csv')
ANI1 = pd.read_csv('ANI1.csv')
ANI1x = pd.read_csv('ANI1x.csv')


def add_data(df_a_e,key): 
    '''
    Add the values of the error for each database to the dataframe with the number of atoms
    by element.

    Parameters
    ----------
    df_a_e : Dataframe with the number of atoms by elements..
    key : Name of the property to be added.

    Returns
    -------
    None.

    '''
    names = [QM9,PC9,ANIE,ANI1, ANI1x]
    Error = []
    for x in names:
        Et = np.array(x[key])
        Error.append(Et)
    names1 = ['QM9','PC9','ANIE','ANI1', 'ANI1x']
    for i in range(0,len(names1)):
        new_key = names1[i] + '_' + key
        df_a_e[new_key] = Error[i]


add_data(df_a_e,'Diff_NN_abi')

df_9_atoms = df_a_e[df_a_e['n_atoms']<=9]
df_9p_atoms = df_a_e[df_a_e['n_atoms']>9]



def calculate_mae(df,atomtype,database):
    '''
    Compute the MAE by element and database. It takes the database created before.

    Parameters
    ----------
    df : Dataframe with errors and number of atoms by element
    atomtype : Name of the element 
    database : Name of the database.

    Returns
    -------
    p : Vector with the number of atoms of a given element and the MAE.

    '''
    unique = np.unique(np.array(df[atomtype]))
    MAE_by_nelement = []
    for i in unique:
        df_x = df[df[atomtype]==i]
        array = np.array(df_x[database])
        mae = np.abs(array.mean())
        MAE_by_nelement.append(mae)
    p = [unique, np.array(MAE_by_nelement)]
    return p

type_of_element = ['n_C', 'n_N', 'n_O']



def get_data(df,key, error):
    '''
    Create an array for each chemical element from the database and compute
    its MAE.

    Parameters
    ----------
    df : Dataframe with the values of the error and the number of atoms of a 
    given chemical element.
    key : Chemical Elemment
    error : Error wish to be studied.

    Returns
    -------
    ad : TYPE
        DESCRIPTION.

    '''
    names1 = ['QM9','PC9','ANIE','ANI1', 'ANI1x']
    ad = []
    for i in range(0,len(names1)):
        new_key = names1[i] + '_' + error
        x = calculate_mae(df,key,new_key)
        ad.append(x)
    return ad

error_by_e_9 = {}
error_by_e_9p = {}
for i in type_of_element:
    error_by_e_9[i] = get_data(df_9_atoms,i,'Diff_NN_abi')
    error_by_e_9p[i] = get_data(df_9p_atoms,i, 'Diff_NN_abi')

#%%
colors = ['firebrick','darkorange','deepskyblue','royalblue','navy']
names1 = ['QM9','PC9','ANIE','ANI1', 'ANI1x']
markers = ['o','o', 'o', 's', 'D']
plt.rc('font', family='sans-serif', size=16)
fig, axs = plt.subplots(nrows=3, ncols=2, figsize=(12,15),sharex='col')

ax1 = axs[0][0]
set1 = error_by_e_9['n_C']
for i in range(0,5):
    ax1.plot(set1[i][0],set1[i][1],"--",marker=markers[i],label=names1[i],color=colors[i])

ax1.set_title(r'$n_{atoms}\leq 9$')
ax1.set_ylabel('MAE(kcal/mol)')
ax1.text(0.7, 0.95, ('C atoms'),transform=ax1.transAxes, fontsize=16,
        verticalalignment='top')
ax1.legend(loc=2)
ax1.yaxis.set_major_formatter(FormatStrFormatter('%10.1f'))
ax2 = axs[0][1]
set1 = error_by_e_9p['n_C']
for i in range(0,5):
    ax2.plot(set1[i][0],set1[i][1],"--",marker=markers[i],label=names1[i],color=colors[i])
ax2.set_title(r'$n_{atoms}>9$')
ax2.legend()
ax3 = axs[1][0]
set1 = error_by_e_9['n_N']
for i in range(0,5):
    ax3.plot(set1[i][0],set1[i][1],"--",marker=markers[i],label=names1[i],color=colors[i])
ax3.text(0.7, 0.95, ('N atoms'),transform=ax3.transAxes, fontsize=16,
        verticalalignment='top')
ax3.yaxis.set_major_formatter(FormatStrFormatter('%10.1f'))
ax3.set_ylabel('MAE(kcal/mol)')
ax4 = axs[1][1]
set1 = error_by_e_9p['n_N']
for i in range(0,5):
        ax4.plot(set1[i][0],set1[i][1],"--",marker=markers[i],label=names1[i],color=colors[i])
ax5 = axs[2][0]
set1 = error_by_e_9['n_O']
for i in range(0,5):
        ax5.plot(set1[i][0],set1[i][1],"--",marker=markers[i],label=names1[i],color=colors[i])
ax5.text(0.7, 0.95, ('O atoms'),transform=ax5.transAxes, fontsize=16,
        verticalalignment='top')
ax5.set_xlabel('Number of atoms')
ax5.set_xlim(left=1)
ax5.yaxis.set_major_formatter(FormatStrFormatter('%10.1f'))
ax5.set_ylabel('MAE(kcal/mol)')
ax6 = axs[2][1]
set1 = error_by_e_9p['n_O']
for i in range(0,5):
        ax6.plot(set1[i][0],set1[i][1],"--",marker=markers[i],label=names1[i],color=colors[i])
ax6.set_xlabel('Number of atoms')
ax6.set_xlim(left=1)   

plt.subplots_adjust(hspace=0)
# plt.savefig('mae_by_elements_tauto.pdf')