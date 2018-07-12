# -*- coding: utf-8 -*-
"""
Created on Tue Jun 26 14:37:47 2018

@author: vlieg
"""
import numpy as np
import functions as func

#%%
def read_logfile(path_logfile):
    path_logfile = func.ChangeExtension(path_logfile,'.log')
    f = open(path_logfile, 'r')
    log = f.readlines()[:]
    f.close()
    
    for n, line in enumerate(log):
        log[n] = line.rstrip()
        
    return log

# Read 16-bit .bin file as 2D image 
def read_bin(path_binfile,slice_nr):
    path_binfile = func.ChangeExtension(path_binfile,'.bin')
    
    path_logfile = func.ChangeExtension(path_binfile,'.log')
    data_log = read_logfile(path_logfile)
    
    # Get 2D dimensions image from log file 
    ypix = np.uint(str.split(data_log[8],' ')[3])
    xpix = np.uint(str.split(data_log[9],' ')[3])
    pix_slice = ypix*xpix
    
    # Open file and read data of 1 slice
    f = open(path_binfile, 'rb')
    
    if slice_nr!=0:
       f.seek(slice_nr*pix_slice*2,0)
       im_slice = np.fromfile(f, count = pix_slice, dtype = '>i2')
    else:
       im_slice = np.fromfile(f, count = pix_slice, dtype = '>i2')  
        
    f.close()
    im_slice = np.reshape(im_slice,[ypix,xpix],'F')

    return im_slice


#%% Read Specified Stack from .bin file 
def get_stack(filepath,stack_nr):
    logfile= read_logfile(filepath)
    
    ypix   = np.uint(str.split(logfile[8],' ')[3])
    xpix   = np.uint(str.split(logfile[9],' ')[3])
    zsteps = logfile[18]
    zsteps = str.split(zsteps," ")[3]
    zsteps = np.uint(str.split(zsteps,",")[0])
    stack=np.zeros([xpix,ypix,zsteps], dtype = np.uint16)
        
    for slice_nr in range(stack_nr*zsteps,stack_nr*zsteps + zsteps):
        stack[:,:,slice_nr-stack_nr*zsteps]=read_bin(filepath,slice_nr)
    
    return stack 
        

