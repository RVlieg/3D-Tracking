# -*- coding: utf-8 -*-
"""
Created on Mon Jul  2 14:52:00 2018

@author: vlieg
"""

import numpy as np
import functions as func
import file_io as file_io
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


#%% Build Stack from 2D slices 
#Get file path 
filepath = 'C:\\Measurement Data\\180618 - PSF GNRs & DOE Pattern\\data_003.bin'


#%% Get Global Peak Coordinates from the entire range of stacks in the file 
threshold = 1.8         #Threshold = threshold*median_stack
mask_size = [11,11,23]   #XYZ in pix 

global_coords = func.get_global_coords(filepath,threshold,mask_size)

#%% Get Local Peak Coordinates by Fitting 3D Gaussian

# Get Measurement Parameters 
logfile = file_io.read_logfile(filepath)
ypix   = np.uint(str.split(logfile[8],' ')[3])
xpix   = np.uint(str.split(logfile[9],' ')[3])
zsteps = logfile[18]
zsteps = str.split(zsteps," ")[3]
zsteps = np.uint(str.split(zsteps,",")[0])
nframes= np.uint(str.split(logfile[7],' ')[2])
nstacks= int(nframes/zsteps) 

#%% Get ROI from the stack using Global Coordinates and Fit 3D Gaussian

# Make 3D Mask
size_x,size_y,size_z=mask_size
mask = func.create_3D_mask(size_x,size_y,size_z,10,10)

# Define Bounds
max_bounds=(65536,65536,size_x,size_y,size_z,size_x,size_y,size_z)
min_bounds=(0,0,0,0,0,0,0,0)
bounds = (min_bounds,max_bounds)
#%%
# Model function to be used to fit to the data:
def gauss_3D(xyz, *p):
    xx,yy,zz = xyz[0],xyz[1],xyz[2]
    A,C,x0,y0,z0,wx,wy,wz = p
    gauss_x = np.square(xx-x0)/(2*np.square(wx))
    gauss_y = np.square(yy-y0)/(2*np.square(wy))
    gauss_z = np.square(zz-z0)/(2*np.square(wz))
    
    return A*np.exp(-1*(gauss_x+gauss_y+gauss_z))+C
        
# Allocate memory for all fit parameters and errors 
num_peaks = len(global_coords[:,0,0])
num_stacks= len(global_coords[0,0,:])

fit_params_all = np.empty([num_peaks,len(max_bounds),num_stacks])
mask_size_all  = np.empty([num_peaks,3,num_stacks])
fit_errors_all = np.empty([num_peaks,len(max_bounds),num_stacks])

#%%
for stack_nr in range(0,len(global_coords[0,0,:])):
        
    # Allocate memory for fit parameters for 1 stack 
    fit_params_stack = np.empty([len(global_coords[:,0,0]),len(max_bounds)])
    fit_errors_stack = np.empty([len(global_coords[:,0,0]),len(max_bounds)])
    mask_size_stack  = np.empty([len(global_coords[:,0,0]),3])
    # Read one stack from .bin file 
    stack = file_io.get_stack(filepath,stack_nr)    
    
    for trace_nr in range(0,len(global_coords[:,0,0])):
        
        # Get ROI
        peak_coords = global_coords[trace_nr,:,stack_nr]
        ROI_stack,mask_size = func.get_ROI_from_stack(filepath,stack,peak_coords,mask)
            
        ##### Fit Data 
        # Get ROI-data from the stack
        data = np.ndarray.flatten(ROI_stack)
        size_x,size_y,size_z = mask_size    
        x = np.arange(0,size_x)
        y = np.arange(0,size_y)
        z = np.arange(0,size_z)
        
        # Create xyz-coordinates 
        xx, yy, zz = np.meshgrid(x, y, z, sparse=False)
        xyz = [np.ndarray.flatten(xx),np.ndarray.flatten(yy), np.ndarray.flatten(zz)]
       
        # Initial Fit Parameters
        in_A  = float(max(data))
        in_C  = float(min(data))
        in_x0 = max(x)/2
        in_y0 = max(y)/2
        in_z0 = max(z)/2
        in_wx = 1.5
        in_wy = 1.5
        in_wz = 3
        
        # p0 is the initial guess for the fitting coefficients
        p0 = [in_A, in_C, in_x0, in_y0, in_z0, in_wx, in_wy, in_wz]
        
        coeff, var_matrix = curve_fit(gauss_3D, xyz, data, p0=p0,bounds=bounds)
        coeff = list(coeff)
        perr = np.sqrt(np.diag(var_matrix))
        
        fit_params_stack[trace_nr,:]=coeff   
        fit_errors_stack[stack_nr,:]=perr
        mask_size_stack[trace_nr,:] = mask_size
        
        
    fit_params_all[:,:,stack_nr]=fit_params_stack 
    fit_errors_all[:,:,stack_nr]=fit_errors_stack
    mask_size_all[:,:,stack_nr] =mask_size_stack 
    
    
#%% Transform Global to Local (sub-pixel) Coordinates 

global_x,global_y,global_z            = fit_params_all[:,2,:],fit_params_all[:,3,:],fit_params_all[:,4,:]
mask_size_x, mask_size_y, mask_size_z = mask_size_all[:,0,:],mask_size_all[:,1,:],mask_size_all[:,2,:]

centered_x = global_x-mask_size_x/2
centered_y = global_y-mask_size_y/2
centered_z = global_z-mask_size_z/2

local_x,local_y,local_z = global_coords[:,0,:]+centered_x, global_coords[:,1,:]+centered_y, global_coords[:,2,:]+centered_z
local_coords = [local_x,local_y,local_z]


