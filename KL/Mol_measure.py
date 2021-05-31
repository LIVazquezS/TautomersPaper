#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 23 15:56:32 2020

@author: L.I.Vazquez-Salazar
@email: litzavazquezs@gmail.com
Class to compute the bond lenght on a molecule from an .xyz file. 

"""


import numpy as np
import itertools

class Mol_measure:
    
    def __init__(self):
        self.radius = []
        self.labels = []
        self.Rtmp = []
        
    def get_indexes(self, x,xs):
        '''
        

        Parameters
        ----------
        x : Label of atom a.
        xs : array of labels of the molecule.

        Returns
        -------
        g_i : Indexes of the atoms in the molecule.

        '''
        g_i = [i for (y, i) in zip(xs, range(len(xs))) if x == y]
        return g_i
    
    def read_molecule(self,file):
        '''
        Read the xyz file format of the molecule.

        Parameters
        ----------
        file : xyz file to be read by the system.

        Returns
        -------
        labels: Array with the labels of the atoms on the molecule
        radius: Array with the radius of the atoms plus a fudge factor for the definition of the bonds
        RTmp: Array of coordinates of each atom.

        '''
        atomic_radius = {'H':0.38, 'C':0.77, 'N':0.75, 'O':0.73, 'F':0.71} # Without fudge factor

        fudge_factor = 0.05

        with open(file) as f:
            contents = f.read().splitlines()
            natom = int(contents[0])
            tag, index, HOMO1, LUMO1, _, _, _, HOMO, LUMO, gap, _, _, Ex, _, _, _, _ = contents[1].split()
                
            self.Rtmp = np.zeros((natom,3))
            self.radius = []
            self.labels = []
            for i in range(natom):
                label, x, y, z, q = contents[2+i].split()
                self.labels.append(label)
                self.Rtmp[i,0] = float(x.replace('*^','e'))
                self.Rtmp[i,1] = float(y.replace('*^','e'))
                self.Rtmp[i,2] = float(z.replace('*^','e'))
                rd = atomic_radius[label] + fudge_factor
                self.radius.append(rd)
    
        return self.labels, self.radius, self.Rtmp
    
    def bond_index(self,a,b):
        '''
        Compute the indexes of the bonds between two atoms 
        and discart possible repetitions. 
        

        Parameters
        ----------
        a : Atom label for atom a.
        b : Atom label for atom b.

        Returns
        -------
        index : Array of indexes to define the bonds.

        '''
        index = []
        x_ind = self.get_indexes(a,self.labels)
        y_ind = self.get_indexes(b,self.labels)
        all_list = [x_ind, y_ind]
        res = list(itertools.product(*all_list))
        for vector in res:
            vector_inverse = vector[::-1]
            if vector_inverse not in index and len(set(vector))==2: 
                index.append(vector)                    
        return index


    def bond_length(self,array):
        '''
        Takes the array that define the bond lenght, checks if it is a valid bond distance
        and return the value of a valid distance.

        Parameters
        ----------
        array : Array of indexes of the molecules

        Returns
        -------
        dist : Bond lenght distance.

        '''
        i = array[0]
        j = array[1]
        r_valid = self.radius[i] + self.radius[j]
        dist = np.linalg.norm(self.Rtmp[i] - self.Rtmp[j])
        if dist > 0.0001 and dist < r_valid:
            return dist
  
    def calculate_bond(self,a,b):
        '''
        Checks for all the combination of the given atoms on the molecule
        compute all the bonds between those atoms on the molecule and return 
        an array with the results. 

        Parameters
        ----------
        a : Atomic label of atom 'A'
        b : Atomic label of atom 'B'

        Returns
        -------
        list_of_bl : Array of bond distance for a pair of atoms on all the molecule.

        '''
        x = self.bond_index(a,b)
        list_of_bl = []
        for element in x:
            bond = self.bond_length(element)
            if bond != None:
                list_of_bl.append(bond)
        return list_of_bl