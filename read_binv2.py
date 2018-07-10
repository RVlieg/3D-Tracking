# -*- coding: utf-8 -*-
"""
Created on Fri Jul  6 15:44:33 2018

@author: Redmar
"""

#%% 
import numpy as np
import functions as func
import file_io as file_io
import bitstring as bitstr
filepath = 'C:\\Measurement Data\\180618 - PSF GNRs & DOE Pattern\\data_003.bin'

#%%

def read_bin(path_binfile,slice_nr):
    path_binfile = func.ChangeExtension(path_binfile,'.bin')
    
    path_logfile = func.ChangeExtension(path_binfile,'.log')
    data_log = file_io.read_logfile(path_logfile)
    
    # Get 2D dimensions image from log file 
    ypix = np.uint(str.split(data_log[8],' ')[3])
    xpix = np.uint(str.split(data_log[9],' ')[3])
    pix_slice = ypix*xpix
    
    # Open file and read data of 1 slice
    f = open(path_binfile, mode='rb')
    
    if slice_nr!=0:
       f.seek(slice_nr*pix_slice*2,0)
       im_slice = np.fromfile(f, dtype = np.uint16, count = pix_slice)
    else:
       im_slice = np.fromfile(f, dtype = np.uint16, count = pix_slice)  
        
    f.close()
    im_slice = np.reshape(im_slice,[ypix,xpix],'F')

    return im_slice
    
aap2 = read_bin(filepath,0)

#%%

path_binfile = func.ChangeExtension(filepath,'.bin')
    
path_logfile = func.ChangeExtension(filepath,'.log')
data_log = file_io.read_logfile(path_logfile)

# Get 2D dimensions image from log file 
ypix = np.uint(str.split(data_log[8],' ')[3])
xpix = np.uint(str.split(data_log[9],' ')[3])
pix_slice = ypix*xpix

# Open file and read data of 1 slice
f = open(path_binfile, 'rb')
slice_nr=1
if slice_nr!=0:
   f.seek(slice_nr*pix_slice*2,0)
   im_slice = np.fromfile(f, count = pix_slice, dtype = '>i2')
else:
   im_slice = np.fromfile(f, count = pix_slice, dtype = np.uint16)  
    
f.close()
im_slice = np.reshape(im_slice,[ypix,xpix],'F')

