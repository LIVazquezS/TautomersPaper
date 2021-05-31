#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  6 14:09:32 2020

@author: L.I. Vazquez-Salazar
@email: litzavazquez@gmail.com
 
File to do a scatter plot of the tautomerization energy. 

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy import stats
import pandas as pd

file_QM9 = pd.read_csv('QM9.csv')
file_PC9 = pd.read_csv('PC9.csv')
file_ANIE= pd.read_csv('ANIE.csv')


# =============================================================================
# 9 atom max

file_QM9_9 = file_QM9[file_QM9['natoms']<=9]
file_PC9_9 = file_PC9[file_PC9['natoms']<=9]
file_ANIE_9 =file_ANIE[file_ANIE['natoms']<=9]

# =============================================================================

file_QM9_9plus = file_QM9[file_QM9['natoms']>9]
file_PC9_9plus = file_PC9[file_PC9['natoms']>9]
file_ANIE_9plus =file_ANIE[file_ANIE['natoms']>9]



def getarray(file):
    '''
    Obtains the pair of values for the tautomerization energy obtained with DFT and 
    the tautomerization energy obtained with NN. 

    Parameters
    ----------
    file : database

    Returns
    -------
    array : Return an array of the pairs of the tautomerization energy.

    '''
    x = file['Diff_abi']
    y = file['Diff_NN']
    array = [x,y]
    return array

# =============================================================================
# 9 atoms
Tauto_QM9_9 = getarray(file_QM9_9)
Tauto_PC9_9 = getarray(file_PC9_9)
Tauto_ANIE_9 = getarray(file_ANIE_9)

# =============================================================================

# =============================================================================
# More than 9
Tauto_QM9_9plus = getarray(file_QM9_9plus)
Tauto_PC9_9plus = getarray(file_PC9_9plus)
Tauto_ANIE_9plus = getarray(file_ANIE_9plus)

# =============================================================================

def stat_quant(array):
    '''
    Compute statistical quantities relevant as MAE, RMSE and Pearson Correlation
    Coefficient ($r^2)
    
    Parameters
    ----------
    array : Pairs of values of the tautomerization energy.

    Returns
    -------
    array_quant : Array of statistical quantities: Pearson correlation coefficient, MAE 
                and RMSE.

    '''
    diff_abi = array[0]
    diff_NN = array[1]
    slope, itercept, r_value, p_value, std_err = stats.linregress(diff_abi, diff_NN)
    r2 = r_value**2
    MSE = np.square(np.subtract(diff_abi,diff_NN)).mean() 
    RMSE = np.sqrt(MSE) 
    MAE = np.abs(np.subtract(diff_abi,diff_NN)).mean()
    array_quant = [r2,MAE,RMSE]
    return array_quant

# =============================================================================
# 9 atoms
QM9_9_stats = stat_quant(Tauto_QM9_9)
PC9_9_stats = stat_quant(Tauto_PC9_9)
ANIE_9_stats = stat_quant(Tauto_ANIE_9)

# =============================================================================

# =============================================================================
# more than 9
QM9_9plus_stats = stat_quant(Tauto_QM9_9plus)
PC9_9plus_stats = stat_quant(Tauto_PC9_9plus)
ANIE_9plus_stats = stat_quant(Tauto_ANIE_9plus)


# =============================================================================
line = np.arange(-150,150)
line1 = np.arange(-150,150)

colors = ['firebrick','darkorange','deepskyblue','royalblue','navy']

fig, axs = plt.subplots(nrows=3, ncols=2, figsize=(10,15),sharex='col',squeeze=False)
plt.rc('font', family='sans-serif', size=16)
mpl.rcParams['axes.linewidth'] = 1.3 
mpl.rcParams['axes.titlepad'] = 10

ax1 = axs[0][0]
# ax1.set_aspect('equal')
ax1.set_title(r'$n_{\rm{atoms}} \leq 9}$ ', size =25)
ax1.plot(Tauto_QM9_9[0],Tauto_QM9_9[1], 'ro',label='B3LYP/6-31G(2df,p)', color=colors[0])
ax1.plot(line,line, color='black')
ax1.text(0.05, 0.97, (r'$\rm{A}$'),transform=ax1.transAxes, fontsize=25,
        verticalalignment='top')
ax1.text(0.5,0.1, (r'$[%.2f,%.2f]$' % (QM9_9_stats[0],QM9_9_stats[1])),transform=ax1.transAxes, fontsize=20,
        verticalalignment='top')
ax1.text(0.4, 0.95, (r'$\rm{QM9}$'),transform=ax1.transAxes, fontsize=20,
        verticalalignment='top')
ax1.set_ylim(-150,150)
ax1.set_ylabel(r'$\Delta E^{\rm{NN}}_{\rm{Tauto}}\rm{(kcal/mol)}$')

ax2 = axs[1][0]
# ax2.set_aspect('equal')
ax2.plot(Tauto_PC9_9[0],Tauto_PC9_9[1], 'ro', label='B3LYP/6-31G(d)',  color=colors[1])
ax2.plot(line,line, color='black')
ax2.text(0.05, 0.97, (r'$\rm{C}$'),transform=ax2.transAxes, fontsize=25,
        verticalalignment='top')
ax2.text(0.65, 0.1, (r'$[%.2f,%.2f]$' % (PC9_9_stats[0],PC9_9_stats[1])),transform=ax2.transAxes, fontsize=14,
         verticalalignment='top')
ax2.text(0.4, 0.95, (r'$\rm{PC9}$'),transform=ax2.transAxes, fontsize=20,
        verticalalignment='top')
ax2.set_ylim(-150,150)
ax2.set_ylabel(r'$\Delta E^{\rm{NN}}_{\rm{Tauto}}\rm{(kcal/mol)}$')

ax3 = axs[2][0]
ax3.plot(Tauto_ANIE_9[0],Tauto_ANIE_9[1], 'ro', label=r'$\omega$B97x/6-31g(d)', color=colors[2])
ax3.plot(line,line, color='black')
ax3.text(0.05, 0.97, (r'$\rm{E}$'),transform=ax3.transAxes, fontsize=25,
        verticalalignment='top')
ax3.text(0.65, 0.1, (r'$[%.2f,%.2f]$' % (ANIE_9_stats[0],ANIE_9_stats[1])),transform=ax3.transAxes, fontsize=14,
        verticalalignment='top')
ax3.text(0.4, 0.95, (r'$\rm{ANI-1E}$'),transform=ax3.transAxes, fontsize=20,
        verticalalignment='top')
ax3.set_xlim(-150,150)
ax3.set_ylim(-150,150)

ax3.set_xlabel(r'$\Delta E^{\rm{DFT}}_{\rm{Tauto}}\rm{(kcal/mol)}$' )
ax3.set_ylabel(r'$\Delta E^{\rm{NN}}_{\rm{Tauto}}\rm{(kcal/mol)}$')


ax6 = axs[0][1]
# ax6.set_aspect('equal')
ax6.set_title(r'$n_{\rm{atoms}} > 9}$ ', size=25)
ax6.plot(Tauto_QM9_9plus[0],Tauto_QM9_9plus[1], 'ro',label='B3LYP/6-31G(2df,p)', color=colors[0])
ax6.plot(line1,line1, color='black')
ax6.text(0.05, 0.97, (r'$\rm{B}$'),transform=ax6.transAxes, fontsize=25,
        verticalalignment='top')
ax6.text(0.65, 0.1, (r'$[%.2f,%.2f]$' % (QM9_9plus_stats[0],QM9_9plus_stats[1])),transform=ax6.transAxes, fontsize=14,
        verticalalignment='top')

ax6.set_ylim(-150,150)
ax6.text(0.4, 0.95, (r'$\rm{QM9}$'),transform=ax6.transAxes, fontsize=20,
        verticalalignment='top')

ax7 = axs[1][1]
# ax7.set_aspect('equal')
ax7.plot(Tauto_PC9_9plus[0],Tauto_PC9_9plus[1], 'ro', label='B3LYP/6-31G(d)', color=colors[1])
ax7.plot(line1,line1, color='black')
ax7.text(0.05, 0.97, (r'$\rm{D}$'),transform=ax7.transAxes, fontsize=25,
        verticalalignment='top')
ax7.text(0.65, 0.1, (r'$[%.2f,%.2f]$' % (PC9_9plus_stats[0],PC9_9plus_stats[1])),transform=ax7.transAxes, fontsize=14,
        verticalalignment='top')
ax7.text(0.4, 0.95, (r'$\rm{PC9}$'),transform=ax7.transAxes, fontsize=20,
        verticalalignment='top')
ax7.set_ylim(-150,150)

ax8 = axs[2][1]
# ax8.set_aspect('equal')
ax8.plot(Tauto_ANIE_9plus[0],Tauto_ANIE_9plus[1], 'ro', label=r'$\omega$B97x/6-31g(d)', color=colors[2])
ax8.plot(line1,line1, color='black')
ax8.text(0.05, 0.97, (r'$\rm{F}$'),transform=ax8.transAxes, fontsize=25,
        verticalalignment='top')
ax8.text(0.65, 0.1, (r'$[%.2f,%.2f]$' % (ANIE_9plus_stats[0],ANIE_9plus_stats[1])),transform=ax8.transAxes, fontsize=14,
        verticalalignment='top')
ax8.text(0.4, 0.95, (r'$\rm{ANI-1E}$'),transform=ax8.transAxes, fontsize=20,
        verticalalignment='top')
ax8.set_xlim(-150,150)
ax8.set_ylim(-150,150)
ax8.set_xlabel(r'$\Delta E^{\rm{DFT}}_{\rm{Tauto}}\rm{(kcal/mol)}$')

plt.subplots_adjust(hspace=0,wspace=0.3)

# plt.tight_layout()
# plt.savefig('Fig1.pdf',bbox_inches='tight')
#plt.show()
