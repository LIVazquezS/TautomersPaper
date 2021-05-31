#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 12:22:25 2021

@author: L.I.Vazquez Salazar
@email: litzavazquezs@gmail.com

Script to plot by number of atoms the MAE on tautomerization energy. 
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition,
                                                  mark_inset)

QM9 = pd.read_csv('QM9_reorder.csv')
PC9 = pd.read_csv('PC9_reorder.csv')
ANIE = pd.read_csv('ANIE_reorder.csv')
ANI1 = pd.read_csv('ANI1_reorder.csv')
ANI1x = pd.read_csv('ANI1x_reorder.csv')

names = [QM9,PC9,ANIE,ANI1, ANI1x]
Errors_Tauto = []
natoms = []
for x in names:
    Et = np.array(x['Diff_NN_abi'])
    nat = np.array(x['n_atoms'])
    Errors_Tauto.append(Et)
    natoms.append(nat)
    
# =============================================================================

def getindex(x,n):
    '''
    Obtaine the indexes of the elements with a given number of atoms 'n'

    Parameters
    ----------
    x : Array with the values of the heavy atoms.
    n : Number of heavy atoms. 

    Returns
    -------
    index : Array of indexes of the atoms with a given number of atoms. 

    '''
    l = len(x)
    index = []
    for i in range(0,l):
        if x[i]==n:
            index.append(i)
    return index

def sep_index(natoms_tauto):
    '''
    Separate the molecules by number of atoms. 

    Parameters
    ----------
    natoms_tauto : Array with the number of atoms.

    Returns
    -------
    v : Vector with the indexes of molecules with more than heavy atoms (first element),
        indexes for molecules with less or equal to nine atoms(second element), number of
        molecules with that number of atoms. 

    '''
    total_natoms = np.unique(natoms_tauto)
    counts = np.unique(natoms_tauto, return_counts=True)
    indexesmax9 = []
    indexes = []
    for y in total_natoms:
            p = getindex(natoms_tauto,y)
            if y <= 9:
                indexesmax9.append(p)
            elif y>=9:
                indexes.append(p)
    v = [indexesmax9, indexes, counts]
    return v


def MAEbyAtoms(Error,index):
    '''
    Computes the MAE for the molecules of a given index. 

    Parameters
    ----------
    Error : Array with all the values of the error.
    index : Array with the indexes of the molecules with a n number of atoms

    Returns
    -------
    errors : Array with the errors by number of atoms. 

    '''
    errors = []
    for y in index: 
        errorx = np.abs(Error[y].mean())
        errors.append(errorx)
    return errors       


# =============================================================================
# Compute indexes and counts

indexes9 = []
indexesmore9 = []
counts = []

for i in range(0,len(natoms)):
    xyz = sep_index(natoms[i])
    x = xyz[0]
    y = xyz[1]
    z = xyz[2]
    indexes9.append(x)
    indexesmore9.append(y)
    counts.append(z)
    
# =============================================================================
# Errors for natoms <= 9
# Error for natoms > 9

Errors_Tauto_9 = []
Errors_Tauto_9plus = []
for i in range(0,len(Errors_Tauto)):
    Errors_Tauto1 = MAEbyAtoms(Errors_Tauto[i],indexes9[i])
    Errors_Tauto2 = MAEbyAtoms(Errors_Tauto[i],indexesmore9[i])
    Errors_Tauto_9.append(Errors_Tauto1)
    Errors_Tauto_9plus.append(Errors_Tauto2)
  
# =============================================================================
#Separate to a max of 9 atoms
total_natoms_max9 = counts[0][0][0:7]
counts_9 = counts[0][1][:7]
total_natoms_more9 = counts[0][0][7:]
counts_9more = counts[0][1][7:]
# =============================================================================

# =============================================================================
# Plotting
names1 = ['QM9','PC9','ANI-1E','ANI-1', 'ANI-1x']
colors = ['firebrick','darkorange','deepskyblue','royalblue','navy']
markers = ['o','o', 'o', 's', 'D']
fig = plt.figure(figsize=(15,5))
# fig.subplots_adjust(right=3)
plt.rc('font', family='sans-serif', size=16)
ax1 = fig.add_subplot(121)

ax1.set_title(r'$n_{\rm{atoms}} \leq 9$ ')
ax1.set_xlabel('Number of atoms')
ax1.bar(total_natoms_max9,counts_9,fill=False)
ax1.text(-0.2, 1.1, ('A'),transform=ax1.transAxes, fontsize=20,
        verticalalignment='top')
ax1.set_ylabel('Frequency')

ax2 = ax1.twinx()
for i in range(0,5):
    ax2.plot(total_natoms_max9,Errors_Tauto_9[i],"--",marker=markers[i],label=names1[i],color=colors[i])
ax2.set_ylim(0,5)
ax2.set_ylabel('MAE(kcal/mol)')
ax2.legend(ncol=6, bbox_to_anchor=(1.2, 1.1),
              loc="lower center", fontsize='small')

ax3 = fig.add_subplot(122)
ax3.set_title(r'$n_{\rm{atoms}} > 9$ ')
ax3.set_xlabel('Number of atoms')
ax3.bar(total_natoms_more9,counts_9more,fill=False)
ax3.text(-0.2, 1.1, ('B'),transform=ax3.transAxes, fontsize=20,
        verticalalignment='top')
ax3.set_ylabel('Frequency')
ax4 = ax3.twinx()
for i in range(0,5):
    ax4.plot(total_natoms_more9,Errors_Tauto_9plus[i],"--",marker=markers[i],label=names1[i],color=colors[i])
ax4.set_ylim(0,50)
ax4.set_ylabel('MAE(kcal/mol)')
ax5 = plt.axes([0,0,1,1])
# Manually set the position and relative size of the inset axes within ax1
ip2 = InsetPosition(ax4,[0.3,0.6,0.4,0.35])
ax5.set_axes_locator(ip2)
for i in range(0,5):
    ax5.plot(total_natoms_more9[0:15],Errors_Tauto_9plus[i][0:15],"--",marker=markers[i],label=names1[i],color=colors[i])
# # Mark the region corresponding to the inset axes on ax1 and draw lines
# # in grey linking the two axes.
ax5.set_ylim(0,7)
ax5.set_xlabel('Number of atoms',fontsize=12)
ax5.set_ylabel('MAE(kcal/mol)',fontsize=12)
mark_inset(ax4, ax5, loc1=3, loc2=4)
plt.subplots_adjust(wspace=0.3)
# plt.tight_layout()
#plt.savefig('maebyatom_onlydft.pdf',bbox_inches='tight')
plt.show()

