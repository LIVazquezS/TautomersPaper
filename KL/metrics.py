#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  4 16:29:08 2020

@author: L.I.Vazquez-Salazar
@email: litzavazquezs@gmail.com

Class to compute the KL divergency for two distributions.
"""

import numpy as np
import pandas as pd
from scipy import stats
from scipy import integrate

class Metrics:
    
    def __init__(self,df,df1):
        self.df = df
        self.df1 = df1
            
    def get_data(self,key,key1):
        '''
        Read .csv file for the reference and target databases

        Parameters
        ----------
        key : Key decribing the values of the bond lenght for the reference 
               database
        key1 : Key describing the values of the bond lenght for the target 
                database

        Returns
        -------
        dp : Values of the bond lenghts for reference and targe.

        '''
        data = self.df[key]
        data_clean = data.dropna()
        data1 = self.df1[key1]
        data1_clean = data1.dropna()
        dp = [data_clean, data1_clean]
        return dp
    
    def support_array(self,array,grid_size=1000):
        '''
        Creates an array for the generation of the Gaussian Kernel distribution 
        from the minimum value of the array to the maximum. The size of the grid
        can be modify.

        Parameters
        ----------
        array : Array of values used for the creation of the kernel distribution.
        grid_size : TYPE, optional. Number of points used for the array.
                    The default is 1000.

        Returns
        -------
        sup_arr : Array of values for the evaluation of the Gaussian Kernel distribution.

        '''
        v_min = np.min(array)-1
        v_max = np.max(array)+1
        sup_arr = np.linspace(v_min,v_max,grid_size)
        return sup_arr

    def support_array_for_int(self,val_min, val_max, grid_size=1000):
        '''
        Creates an array for the generation of the Gaussian Kernel distribution 
        from the minimum value of the array to the maximum. The size of the grid
        can be modify. NOTE: THis function is used for integration. 

        Parameters
        ----------
        array : Array of values used for the creation of the kernel distribution.
        grid_size : TYPE, optional. Number of points used for the array.
                    The default is 1000.

        Returns
        -------
        sup_arr : Array of values for the evaluation of the Gaussian Kernel distribution.

        '''
        v_min = val_min 
        v_max = val_max
        sup_arr1 = np.linspace(v_min,v_max,grid_size)
        return sup_arr1     
    
    def kdg(self,array):
        '''
        Gaussian Kernel Distribution

        Parameters
        ----------
        array : Values of the bond distances used for the creation of the 
                Gaussian Kernel distribution.

        Returns
        -------
        kernel_val : Values of the Gaussian Kernel distribution and the values
                    evaluated.

        '''
        kernel = stats.gaussian_kde(array)
        x_eval = self.support_array(array) 
        kernel_val = [x_eval,kernel(x_eval)]
        return kernel_val 
    
    def kdg_int(self,array,val_min,val_max):
        '''
        Gaussian Kernel Distribution for integration.

        Parameters
        ----------
        array : Values of the bond distances used for the creation of the 
                Gaussian Kernel distribution.

        Returns
        -------
        kernel_val : Values of the Gaussian Kernel distribution and the values
                    evaluated.

        '''
        kernel = stats.gaussian_kde(array)
        x_eval = self.support_array_for_int(val_min,val_max)
        kernel_val1 = [x_eval, kernel(x_eval)]
        return kernel_val1 
     

    def prep_KL(self,dp):
        '''
        Take the values of the bond lenght distributions and add a 1 to 
        normalize the function and avoid discontinuities on the evaluation
        of the KL divergency. Additionally, check that the number of points
        of the Kernel distribution is the same to assure a correct calculation. 
        Then it computes the values of the reference distribution multiplied by 
        the values of the logarithmic coefficient between reference and target distribution.

        Parameters
        ----------
        dp : Values of the bond lenghts.

        Returns
        -------
        x_val : Values evaluated on the Gaussian Kernel Distribution
        prod_to_int : Values to integrate.

        '''
        p = self.kdg(dp[0])
        q = self.kdg(dp[1])
        p_norm = p[1] + 1
        q_norm = q[1] + 1
        x_val = q[0]
        if len(p_norm) == len(q_norm):
            log_coeff = np.log(p_norm) - np.log(q_norm)
            prod_to_int = p_norm*log_coeff  
            return x_val, prod_to_int 
    
    def KL_divergence_cum(self,dp):
        '''
        Implementation of Kullback-Leibler divergence as a cummulative integral. 
        It obtains the value on all the points of x.
        IMPORTANT NOTE: the divergence is not symmetric so it is not the same to take D_{kl}(P|Q)
        than D_{kl}(Q|P) 

        Parameters
        ----------
        dp : Values of the bond lenghts

        Returns
        -------
        x_val : Points evaluated on the KL divergence.
        KL : Value of the KL divergency.

        '''

        x_val, values = self.prep_KL(dp)
        KL = integrate.cumtrapz(values, x=x_val)
        return x_val, KL
    
    def KL_divergence_all(self,dp):
        '''
        Implementation of Kullback-Leibler divergence as a integral in all the space.
        It returns a scalar value that is the overall value of the integral.
        IMPORTANT NOTE: the divergence is not symmetric so it is not the same to take D_{kl}(P|Q)
        than D_{kl}(Q|P) 

        Parameters
        ----------
        dp : Values of the bond lenghts for target and reference databases

        Returns
        -------
        KL : Value of the KL divergency.

        '''
        x_val, values = self.prep_KL(dp)
        KL = integrate.trapz(values,x=x_val)
        return KL
                            
    def KL_divergence_int(self, dp,val_min,val_max):
        '''
        Implementation of Kullback-Leibler divergence as a integral on a define interval
        It returns a scalar value that is the value of KL on that section.
        IMPORTANT NOTE: the divergence is not symmetric so it is not the same to take D_{kl}(P|Q)
        than D_{kl}(Q|P) 

        Parameters
        ----------
        dp : Values of the bond lenghts for target and reference databases.
        val_min : Minimum value to consider for the integration.
        val_max : Maximum value to consider for the integration.

        Returns
        -------
        KL_int : Value of the KL divergency.

        '''
        
        p = self.kdg_int(dp[0],val_min,val_max)
        q = self.kdg_int(dp[1],val_min,val_max)
        p_norm = p[1] + 1
        q_norm = q[1] + 1
        x_val = q[0]
        if len(p_norm) == len(q_norm):
            log_coeff = np.log(p_norm) - np.log(q_norm)
            prod_to_int = p_norm*log_coeff
        KL_int = integrate.trapz(prod_to_int,x_val)
        return KL_int
        
                            

        
        