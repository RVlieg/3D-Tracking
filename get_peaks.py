# -*- coding: utf-8 -*-
"""
Created on Mon Jul  2 14:52:00 2018

@author: vlieg
"""

import numpy as np
import functions as func
import file_io as file_io
import matplotlib.pyplot as plt


#%% Build Stack from 2D slices 

#Get file path 
filepath = 'C:\\Measurement Data\\180618 - PSF GNRs & DOE Pattern\\data_003.bin'


#%% Get Global Peak Coordinates from the entire range of stacks in the file 

threshold = 1.8         #Threshold = threshold*median_stack
mask_size = [11,11,23]   #XYZ in pix 

global_coords = func.get_global_coords(filepath,threshold,mask_size)



#%% Get Local Peak Coordinates by Fitting 3D Gaussian

logfile = file_io.read_logfile(filepath)
ypix   = np.uint(str.split(logfile[8],' ')[3])
xpix   = np.uint(str.split(logfile[9],' ')[3])
zsteps = logfile[18]
zsteps = str.split(zsteps," ")[3]
zsteps = np.uint(str.split(zsteps,",")[0])
nframes= np.uint(str.split(logfile[7],' ')[2])
nstacks= int(nframes/zsteps) 

# Read stack from .bin file 
stack = file_io.get_stack(filepath,0)

