#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  9 17:23:58 2020

@author: L.I.VazquezSalazar
@email: litzavazquezs@gmail.com

Script to plot a Gaussian Kernel distribution and violin plots for the error
on the prediction of the tautomerization energy. 
"""

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib
from matplotlib.ticker import FormatStrFormatter

file_QM9 = pd.read_csv('QM9.csv')
file_PC9 = pd.read_csv('PC9.csv')
file_ANIE = pd.read_csv('ANIE.csv')
file_ANI1 = pd.read_csv('ANI1.csv')
file_ANI1x = pd.read_csv('ANI1x.csv')

file_QM9_9 = file_QM9[file_QM9['natoms']<=9]
file_PC9_9 = file_PC9[file_PC9['natoms']<=9]
file_ANIE_9 = file_ANIE[file_ANIE['natoms']<=9]
file_ANI1_9 = file_ANI1[file_ANI1['natoms']<=9]
file_ANI1x_9 = file_ANI1x[file_ANI1x['natoms']<=9]

file_QM9_9plus = file_QM9[file_QM9['natoms']>9]
file_PC9_9plus = file_PC9[file_PC9['natoms']>9]
file_ANIE_9plus = file_ANIE[file_ANIE['natoms']>9]
file_ANI1_9plus = file_ANI1[file_ANI1['natoms']>9]
file_ANI1x_9plus = file_ANI1x[file_ANI1x['natoms']>9]

def filt(array):
    '''
    Filter the values to the 95 percentile of the distribution to show a more
    compact distribution.

    Parameters
    ----------
    array : Values of the error

    Returns
    -------
    new_a : Array with the values filtered to the 95% percentil.

    '''
    mean = np.percentile(array,95)
    new_array =[]
    for i in array:
        if i <= mean:
            new_array.append(i)
        else:
            new_array.append(0)
    new_a = np.array(new_array)
    return new_a


def getdf9(key):
    '''
    Creates a dictionary for a property on the datasets using the filter function. This is for the molecules with
    less or equal to 9 atoms. 

    Parameters
    ----------
    key : Value of the property should be a string.

    Returns
    -------
    df : dataframe with the values of a property for all the databases.

    '''
    E_ANI1x = filt(np.abs(file_ANI1x_9[key])) 
    E_ANI1 = filt(np.abs(file_ANI1_9[key]))
    E_ANIE = filt(np.abs(file_ANIE_9[key]))
    E_PC9 = filt(np.abs(file_PC9_9[key]))
    E_QM9 = filt(np.abs(file_QM9_9[key]))

    dictionary = {'QM9':E_QM9, 'PC9':E_PC9, 'ANI-1E':E_ANIE, 'ANI-1':E_ANI1, 'ANI-1x':E_ANI1x }
    df = pd.DataFrame(dictionary)
    return df


def getdf9plus(key):
    '''
    Creates a dictionary for a property on the datasets using the filter function.
    This is for the molecules with more than 9 atoms.

    Parameters
    ----------
    key : Value of the property should be a string.

    Returns
    -------
    df : dataframe with the values of a property for all the databases.

    '''
    E_ANI1x = filt(np.abs(file_ANI1x_9plus[key])) 
    E_ANI1 = filt(np.abs(file_ANI1_9plus[key]))
    E_ANIE = filt(np.abs(file_ANIE_9plus[key]))
    E_PC9 = filt(np.abs(file_PC9_9plus[key]))
    E_QM9 = filt(np.abs(file_QM9_9plus[key]))
    dictionary = {'QM9':E_QM9, 'PC9':E_PC9, 'ANI-1E':E_ANIE, 'ANI-1':E_ANI1, 'ANI-1x':E_ANI1x }
    df = pd.DataFrame(dictionary)
    return df

def getdf91(key):
   '''
    Creates a dictionary for a property on the datasets without the filter function.
    This is for the molecules with more than 9 atoms.

    Parameters
    ----------
    key : Value of the property should be a string.

    Returns
    -------
    df : dataframe with the values of a property for all the databases.

   '''
    E_ANI1x = file_ANI1x_9[key]
    E_ANI1 = file_ANI1_9[key]
    E_ANIE = file_ANIE_9[key]
    E_PC9 = file_PC9_9[key]
    E_QM9 = file_QM9_9[key]

    dictionary = {'QM9':E_QM9, 'PC9':E_PC9, 'ANI-1E':E_ANIE, 'ANI-1':E_ANI1, 'ANI-1x':E_ANI1x }
    df = pd.DataFrame(dictionary)
    return df


def getdf9plus1(key):
    '''
    Creates a dictionary for a property on the datasets without the filter function.
    This is for the molecules with more than 9 atoms.

    Parameters
    ----------
    key : Value of the property should be a string.

    Returns
    -------
    df : dataframe with the values of a property for all the databases.

    '''
    E_ANI1x = file_ANI1x_9plus[key]
    E_ANI1 = file_ANI1_9plus[key]
    E_ANIE = file_ANIE_9plus[key]
    E_PC9 = file_PC9_9plus[key]
    E_QM9 = file_QM9_9plus[key]

    dictionary = {'QM9':E_QM9, 'PC9':E_PC9, 'ANI-1E':E_ANIE, 'ANI-1':E_ANI1, 'ANI-1x':E_ANI1x }
    df = pd.DataFrame(dictionary)
    return df

Tauto_energy9 = getdf9('Diff_NN_abi')
Tauto_energy9plus = getdf9plus('Diff_NN_abi')

Tauto_energy91 = getdf91('Diff_NN_abi')
Tauto_energy91plus = getdf9plus1('Diff_NN_abi')

colors = ['firebrick','darkorange','deepskyblue','royalblue','navy']
colors_rgba = []
for i in colors:
    x = matplotlib.colors.to_rgba(i)
    colors_rgba.append(x)

fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(11,10),sharex=False)
plt.rc('font', family='sans-serif', size=16)
ax1 = sns.kdeplot(data=Tauto_energy91, ax=axs[0][0], palette=colors_rgba,common_norm=False )
ax1.set_xlim(-20,20)
ax1.set_ylim(0,0.6)
ax1.text(-0.3, 1.1, ('A'),transform=ax1.transAxes, fontsize=25,
        verticalalignment='top')
ax1.yaxis.set_major_formatter(FormatStrFormatter('%10.2f'))
ax1.set_xlabel(r'$\Delta E^{\rm{DFT}}_{\rm{Tauto}}-\Delta E^{\rm{NN}}_{\rm{Tauto}}\rm{(kcal/mol)}$')
ax1.set_title(r'$n_{\rm{atoms}} \leq 9}$ ')

ax2 = sns.kdeplot(data=Tauto_energy91plus, ax=axs[0][1], palette=colors_rgba,common_norm=False)
ax2.text(-0.3, 1.1, ('B'),transform=ax2.transAxes, fontsize=25,
        verticalalignment='top')
ax2.set_xlim(-20,20)
ax2.set_ylim(0,0.6)
ax2.yaxis.set_major_formatter(FormatStrFormatter('%10.2f'))
ax2.set_xlabel(r'$\Delta E^{\rm{DFT}}_{\rm{Tauto}}-\Delta E^{\rm{NN}}_{\rm{Tauto}}\rm{(kcal/mol)}$')
ax2.set_title(r'$n_{\rm{atoms}} > 9$ ')

ax3 = sns.violinplot(data=Tauto_energy9, ax=axs[1][0], inner="box", cut=0, palette=colors_rgba )
sns.boxplot(data=Tauto_energy9, ax=axs[1][0], whis=[0.5,0.95],showbox=False, showmeans=False,showfliers=False)
ax3.set_ylim(0,20)
ax3.text(-0.3, 1.1, ('C'),transform=ax3.transAxes, fontsize=25,
        verticalalignment='top')
ax3.set_ylabel(r'$\Delta E^{\rm{DFT}}_{\rm{Tauto}}-\Delta E^{\rm{NN}}_{\rm{Tauto}}\rm{(kcal/mol)}$')

ax4 = sns.violinplot(data=Tauto_energy9plus, ax=axs[1][1],cut=0, palette=colors_rgba)
sns.boxplot(data=Tauto_energy9plus, ax=axs[1][1], whis=[0.5,0.95],showbox=False, showmeans=False,showfliers=False)
ax4.text(-0.3, 1.1, ('D'),transform=ax4.transAxes, fontsize=25,
        verticalalignment='top')
ax4.set_ylim(0,20)
ax4.set_ylabel(r'$\Delta E^{\rm{DFT}}_{\rm{Tauto}}-\Delta E^{\rm{NN}}_{\rm{Tauto}}\rm{(kcal/mol)}$')
plt.tight_layout() 
#plt.savefig('mixed_onlyDFT.pdf',bbox_inches='tight')
#plt.show()