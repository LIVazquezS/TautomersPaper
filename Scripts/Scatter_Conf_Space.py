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
from scipy import stats
import pandas as pd
import matplotlib as mpl

file_ANI1= pd.read_csv('ANI.csv')
file_ANI1x = pd.read_csv('ANI1x.csv')


file_ANI1_9 = file_ANI1[file_ANI1['natoms']<=9]
file_ANI1x_9 = file_ANI1x[file_ANI1x['natoms']<=9]



file_ANI1_9plus = file_ANI1[file_ANI1['natoms']>9]
file_ANI1x_9plus = file_ANI1x[file_ANI1x['natoms']>9]



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

Tauto_ANI1_9 = getarray(file_ANI1_9)
Tauto_ANI1x_9 = getarray(file_ANI1x_9)


Tauto_ANI1_9plus = getarray(file_ANI1_9plus)
Tauto_ANI1x_9plus = getarray(file_ANI1x_9plus)



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
ANI1_9_stats = stat_quant(Tauto_ANI1_9)
ANI1x_9_stats = stat_quant(Tauto_ANI1x_9)


ANI1_9plus_stats = stat_quant(Tauto_ANI1_9plus)
ANI1x_9plus_stats = stat_quant(Tauto_ANI1x_9plus)


line = np.arange(-100,60)
line1 = np.arange(-150,150)
colors = ['firebrick','darkorange','deepskyblue','royalblue','navy']


fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(10,10),sharex='col',squeeze=False)
plt.rc('font', family='sans-serif', size=16)
mpl.rcParams['axes.linewidth'] = 1.3 
mpl.rcParams['axes.titlepad'] = 10
ax4 = axs[0][0]
ax4.set_title(r'$n_{\rm{atoms}} \leq 9}$ ', size=25)
# ax4.set_aspect('equal')
ax4.plot(Tauto_ANI1_9[0],Tauto_ANI1_9[1], 'ro',label=r'$\omega$B97x/6-31g(d)', color=colors[3])
ax4.plot(line,line, color='black')
ax4.text(0.05, 0.97, (r'$\rm{A}$'),transform=ax4.transAxes, fontsize=25,
        verticalalignment='top')
ax4.text(0.65, 0.1, (r'$[%.2f,%.2f]$' % (ANI1_9_stats[0],ANI1_9_stats[1])),transform=ax4.transAxes, fontsize=14,
        verticalalignment='top')
ax4.text(0.3, 0.95, (r'$\rm{ANI-1}$'),transform=ax4.transAxes, fontsize=20,
        verticalalignment='top')
ax4.set_ylim(-60,60)
ax4.set_ylabel(r'$\Delta E^{\rm{NN}}_{\rm{Tauto}}\rm{(kcal/mol)}$')

ax5 = axs[1][0]
# ax5.set_aspect('equal')
ax5.plot(Tauto_ANI1x_9[0],Tauto_ANI1x_9[1], 'ro', label=r'$\omega$B97x/6-31g(d)', color=colors[4])
ax5.plot(line,line, color='black')
ax5.text(0.05, 0.95, (r'$\rm{C}$'),transform=ax5.transAxes, fontsize=25,
        verticalalignment='top')
ax5.text(0.65, 0.1, (r'$[%.2f,%.2f]$' % (ANI1x_9_stats[0],ANI1x_9_stats[1])),transform=ax5.transAxes, fontsize=14,
        verticalalignment='top')
ax5.text(0.3, 0.95, (r'$\rm{ANI-1x}$'),transform=ax5.transAxes, fontsize=20,
        verticalalignment='top')
ax5.set_xlim(-60,60)
ax5.set_ylim(-60,60)
ax5.set_xlabel(r'$\Delta E^{\rm{DFT}}_{\rm{Tauto}}\rm{(kcal/mol)}$')
ax5.set_ylabel(r'$\Delta E^{\rm{NN}}_{\rm{Tauto}}\rm{(kcal/mol)}$')

ax9 = axs[0][1]
# ax9.set_aspect('equal')
ax9.set_title(r'$n_{\rm{atoms}} > 9}$ ', size=25)
ax9.plot(Tauto_ANI1_9plus[0],Tauto_ANI1_9plus[1], 'ro',label=r'$\omega$B97x/6-31g(d)', color=colors[3])
ax9.plot(line1,line1, color='black')
ax9.text(0.05, 0.97, (r'$\rm{B}$'),transform=ax9.transAxes, fontsize=25,
        verticalalignment='top')
ax9.text(0.65, 0.1, (r'$[%.2f,%.2f]$' % (ANI1_9plus_stats[0],ANI1_9plus_stats[1])),transform=ax9.transAxes, fontsize=14,
        verticalalignment='top')
ax9.text(0.3, 0.95, (r'$\rm{ANI-1}$'),transform=ax9.transAxes, fontsize=20,
        verticalalignment='top')
ax9.set_ylim(-150,150)

ax10 = axs[1][1]
# ax10.set_aspect('equal')
ax10.plot(Tauto_ANI1x_9plus[0],Tauto_ANI1x_9plus[1], 'ro', label=r'$\omega$B97x/6-31g(d)', color=colors[4])
ax10.plot(line1,line1, color='black')
ax10.text(0.05, 0.95, (r'$\rm{D}$'),transform=ax10.transAxes, fontsize=25,
        verticalalignment='top')
ax10.text(0.65, 0.1, (r'$[%.2f,%.2f]$' % (ANI1x_9plus_stats[0],ANI1x_9plus_stats[1])),transform=ax10.transAxes, fontsize=14,
        verticalalignment='top')
ax10.text(0.3, 0.95, (r'$\rm{ANI-1x}$'),transform=ax10.transAxes, fontsize=20,
        verticalalignment='top')
ax10.set_xlim(-150,150)
ax10.set_ylim(-150,150)
ax10.set_xlabel(r'$\Delta E^{\rm{DFT}}_{\rm{Tauto}}\rm{(kcal/mol)}$')

plt.subplots_adjust(hspace=0,wspace=0.3)
# plt.tight_layout()
#plt.savefig('Scatter_confspace_onlydft.pdf',bbox_inches='tight')
#plt.show()

