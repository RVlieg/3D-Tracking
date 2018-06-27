# -*- coding: utf-8 -*-
"""
Created on Tue Jun 26 14:37:47 2018

@author: vlieg
"""
import numpy as np

def read_logfile(path_logfile):

    f = open(path_logfile, 'r')
    log = f.readlines()[:]
    f.close()
    
    for n, line in enumerate(log):
        log[n] = line.rstrip()
        
    return log



# Read 16-bit .bin file as 2D image 
path_binfile = 'C:\\Measurement Data\\180618 - PSF GNRs & DOE Pattern\\data_003.bin'


data = np.fromfile(path_binfile, dtype=np.uint16) 
read_logfile()


2D_slice = np.reshape 

