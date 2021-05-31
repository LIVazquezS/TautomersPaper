#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 17:31:29 2021

@author: L.I.Vazquez-Salazar
@email: litzavazquezs@gmail.com

Script to analyze the error on the prediction of the energy of a molecule by 
number of bonds on a molecule. It need two extra file with the number of bonds
by type on molecules of type A and B. 
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.ticker import FormatStrFormatter

#This can be improved...To do later.
df_bonds_A = pd.read_csv('number_of_bonds_molA.csv')
df_bonds_B = pd.read_csv('number_of_bonds_molB.csv')
E_ANIE = pd.read_csv('ANIE.csv')
E_PC9 = pd.read_csv('PC9.csv')
E_QM9 = pd.read_csv('QM9.csv')
E_ANI1 = pd.read_csv('ANI1.csv')
E_ANI1x = pd.read_csv('ANI1x.csv')

df_bonds_A['ANIE'] = E_ANIE['Diff_MolA']
df_bonds_A['PC9']= E_PC9['Diff_MolA']
df_bonds_A['QM9']= E_QM9['Diff_MolA']
df_bonds_A['ANI1'] = E_ANI1['Diff_MolA']
df_bonds_A['ANI1x'] = E_ANI1x['Diff_MolA']


df_bonds_B['ANIE'] = E_ANIE['Diff_molB']
df_bonds_B['PC9'] = E_PC9['Diff_molB']
df_bonds_B['QM9'] = E_QM9['Diff_molB']
df_bonds_B['ANI1'] = E_ANI1['Diff_molB']
df_bonds_B['ANI1x'] = E_ANI1x['Diff_molB']


df_9_A = df_bonds_A[df_bonds_A['natoms']<=9]
df_9_B = df_bonds_B[df_bonds_B['natoms']<=9]

df_9plus_A = df_bonds_A[df_bonds_A['natoms']>9]
df_9plus_B = df_bonds_B[df_bonds_B['natoms']>9]


df_9 = pd.concat([df_9_A,df_9_B])
df_9p = pd.concat([df_9plus_A,df_9plus_B])

def calculate_mae(df,bondtype,database):
    '''
    Compute the MAE by bond and database. It takes the database created before.

    Parameters
    ----------
    df : Dataframe with errors and number of bonds 
    bondtype : Name of bond
    database : Name of the database.

    Returns
    -------
    p : Vector with the number of bonds of a given bond and the MAE.

    '''
    unique = np.unique(np.array(df[bondtype]))
    MAE_by_nbonds = []
    for i in unique:
        df_x = df[df[bondtype]==i]
        array = np.array(df_x[database])
        mae = np.abs(array.mean())
        MAE_by_nbonds.append(mae)
    p = [unique, np.array(MAE_by_nbonds)]
    return p

types_of_bonds = ['CH', 'CN', 'CO','CC','NN','NO','NH','OH']


def get_data(df,key):
    '''
    Compute the MAE for a given bond and a given dataframe

    Parameters
    ----------
    df : Dataframe. Contains all the values of the error and the number of bonds
    key : String. Identifier of the bond to be described

    Returns
    -------
    ad : Array with errors for all the databases.

    '''
    x = calculate_mae(df,key,'QM9')
    y = calculate_mae(df,key,'PC9')
    z = calculate_mae(df,key,'ANIE')
    x1 = calculate_mae(df,key,'ANI1')
    x2 = calculate_mae(df, key, 'ANI1x')
    ad = [x,y,z,x1,x2]
    return ad

bonds_9 = {}
bonds_9p = {}
for i in types_of_bonds:
    bonds_9[i] = get_data(df_9,i)
    bonds_9p[i] = get_data(df_9p,i)

def get_num_of_mols(df):
    '''
    Counts the number of molecules with a given number of bonds.

    Parameters
    ----------
    df : Data of number of bonds and error.

    Returns
    -------
    data : Repetitions for each type of bond.

    '''
    data = []
    for i in types_of_bonds:
        unique = np.unique(np.array(df[i]),return_counts=True)
        data.append(unique)
    return data

mol_num9 = get_num_of_mols(df_9)
mol_num9p = get_num_of_mols(df_9p)


#%%
colors = ['firebrick','darkorange','deepskyblue','royalblue','navy']
names1 = ['QM9','PC9','ANI-1E','ANI-1', 'ANI-1x']
markers = ['o','o', 'o', 's', 'D']
plt.rc('font', family='sans-serif', size=16)
fig, axs = plt.subplots(nrows=3, ncols=2, figsize=(15,15),sharex='col')

ax1 = axs[0][0]
set1 = bonds_9['CC']
for i in range(0,5):
    ax1.plot(set1[i][0],set1[i][1],"--",marker=markers[i],label=names1[i],color=colors[i])

ax1.set_title(r'$n_{\rm{atoms}}\leq 9$')
ax1.set_ylabel('MAE(kcal/mol)')
ax1.text(0.7, 0.95, ('C-C bonds'),transform=ax1.transAxes, fontsize=16,
        verticalalignment='top')
ax1.text(0.05, 0.95, ('A'),transform=ax1.transAxes, fontsize=24,
        verticalalignment='top')
ax1.legend(ncol=5, bbox_to_anchor=(1.2, 1.1),
              loc="lower center", fontsize='small')
ax1.set_ylim(bottom=0)
ax1bis = ax1.twinx()
ax1bis.bar(mol_num9[3][0],mol_num9[3][1],fill=False)
ax1bis.set_ylim(bottom=0)
ax2 = axs[0][1]
set1 = bonds_9p['CC']
for i in range(0,5):
    ax2.plot(set1[i][0],set1[i][1],"--",marker=markers[i],label=names1[i],color=colors[i])
ax2.set_title(r'$n_{\rm{atoms}}>9$')
ax2.text(0.05, 0.95, ('B'),transform=ax2.transAxes, fontsize=24,
        verticalalignment='top')
ax2.set_ylim(bottom=0)
ax2bis = ax2.twinx()
ax2bis.bar(mol_num9p[3][0],mol_num9p[3][1],fill=False)
ax2bis.set_ylabel('Frequency')


ax3 = axs[1][0]
set1 = bonds_9['CO']
for i in range(0,5):
    ax3.plot(set1[i][0],set1[i][1],"--",marker=markers[i],label=names1[i],color=colors[i])
ax3.text(0.7, 0.95, ('C-O bonds'),transform=ax3.transAxes, fontsize=16,
        verticalalignment='top')
ax3.text(0.05, 0.95, ('C'),transform=ax3.transAxes, fontsize=24,
        verticalalignment='top')
ax3.set_ylabel('MAE(kcal/mol)')
ax3.set_ylim(bottom=0)
ax3bis = ax3.twinx()
ax3bis.bar(mol_num9[2][0],mol_num9[2][1],fill=False)


ax4 = axs[1][1]
set1 = bonds_9p['CO']
for i in range(0,5):
        ax4.plot(set1[i][0],set1[i][1],"--",marker=markers[i],label=names1[i],color=colors[i])
 
ax4.text(0.05, 0.95, ('D'),transform=ax4.transAxes, fontsize=24,
        verticalalignment='top')
ax4.set_ylim(bottom=0)        
ax4bis = ax4.twinx()
ax4bis.bar(mol_num9p[2][0],mol_num9p[2][1],fill=False)
ax4bis.set_ylabel('Frequency')        

ax5 = axs[2][0]
set1 = bonds_9['CN']
for i in range(0,5):
        ax5.plot(set1[i][0],set1[i][1],"--",marker=markers[i],label=names1[i],color=colors[i])
ax5.text(0.7, 0.95, ('C-N bonds'),transform=ax5.transAxes, fontsize=16,
        verticalalignment='top')
ax5.text(0.05, 0.95, ('E'),transform=ax5.transAxes, fontsize=24,
        verticalalignment='top')
ax5.set_xlabel('Number of bonds')
ax5.set_xlim(left=0)
ax5.yaxis.set_major_formatter(FormatStrFormatter('%i'))
ax5.set_ylabel('MAE(kcal/mol)')
ax5bis = ax5.twinx()
ax5bis.bar(mol_num9[1][0],mol_num9[1][1],fill=False)


ax6 = axs[2][1]
set1 = bonds_9p['CN']
for i in range(0,5):
        ax6.plot(set1[i][0],set1[i][1],"--",marker=markers[i],label=names1[i],color=colors[i])
ax6.set_xlabel('Number of bonds')
ax6.text(0.1, 0.95, ('F'),transform=ax6.transAxes, fontsize=24,
        verticalalignment='top')
ax6.set_xlim(left=0)   
ax6bis = ax6.twinx()
ax6bis.bar(mol_num9p[1][0],mol_num9p[1][1],fill=False)
ax6bis.set_ylabel('Frequency')
plt.subplots_adjust(hspace=0,wspace=0.3)
#plt.savefig('num_of_bonds_mae_C.pdf',bbox_inches='tight')
# plt.tight_layout()
