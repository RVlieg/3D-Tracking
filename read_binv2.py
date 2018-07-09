# -*- coding: utf-8 -*-
"""
Created on Fri Jul  6 15:44:33 2018

@author: Redmar
"""

#%% 
import numpy as np
filepath = 'C:\\Users\\Redmar\\Desktop\\data_003.bin'


with open(filepath, "rb") as binary_file:
    aap = binary_file.read()
    
    aap = np.array(aap)
    binary_file.close