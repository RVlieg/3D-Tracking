# -*- coding: utf-8 -*-
"""
Created on Thu Jun 28 14:57:53 2018

@author: vlieg
"""
import numpy as np
import tkinter as tk
import file_io as file_io
from scipy.optimize import curve_fit
from tkinter import filedialog

#%% Change Extension of a String
def ChangeExtension(path_file, extension):
    path_old = path_file    
    extension_index = str.find(path_old, '.')    
    if str.find(extension,'.') != 0:
        extension = '.' + extension        
    path_new = path_old[0:extension_index] + extension

    return path_new
    

#%% Get filepath from a pop-up window
def get_path():
    root = tk.Tk()
    root.withdraw()

    filename = filedialog.askopenfilename()
    filename = filename.replace('/','\\')

    return filename 


#%% Get max value in array with corresponding indices
def get_max(array):
    indices = np.unravel_index(np.argmax(array), np.shape(array))
    max_intensity = array[indices]
    indices = np.int16(indices)    
    return max_intensity, indices
    
    
#%% Get min value in array with corresponding indices
def get_min(array):
    indices = np.unravel_index(np.argmin(array), np.shape(array))
    max_intensity = array[indices]
    indices = np.int16(indices)
    return max_intensity, indices
    
    
#%% Get Peaks from Stack 
def create_3D_mask(size_x,size_y,size_z,amp,sigma):
    # Create 3D Gaussian mask 
    
    A = amp   # amplitude
    S = sigma # sigma
    
    x = np.arange(0,size_x,1)
    y = np.arange(0,size_y,1)
    z = np.arange(0,size_z,1)
    mid_indx = np.uint16(np.floor(len(x)/2))
    mid_indy = np.uint16(np.floor(len(y)/2))
    mid_indz = np.uint16(np.floor(len(z)/2))
    
    # Create xyz-coordinates 
    xx, yy, zz = np.meshgrid(x, y, z, sparse=False)
    
    # Calculate 3D Gaussian
    gauss_x = np.square(xx-x[mid_indx])/(2*np.square(S))
    gauss_y = np.square(yy-y[mid_indy])/(2*np.square(S))
    gauss_z = np.square(zz-z[mid_indz])/(2*np.square(S))
    gauss_xyz = A*np.exp(-1*(gauss_x+gauss_y+gauss_z))
    
    # Get Mask in boolean
    mask_xyz = gauss_xyz<np.median(gauss_xyz)
    
    return mask_xyz

#%% Get ROI from stack 

def get_ROI_from_stack(filepath,stack,peak_coords,mask):
    
    # Get file parameters     
    logfile = file_io.read_logfile(filepath)

    ypix = np.uint(str.split(logfile[8],'=')[1])
    xpix = np.uint(str.split(logfile[9],'=')[1])
    zsteps = logfile[18]
    zsteps = str.split(zsteps,"=")[1]
    zsteps = np.uint(str.split(zsteps,",")[0])   
    size_x,size_y,size_z=mask.shape
    
    x_range = np.int16([peak_coords[0]-np.floor(size_x/2),peak_coords[0]+np.floor(size_x/2)])
    y_range = np.int16([peak_coords[1]-np.floor(size_y/2),peak_coords[1]+np.floor(size_y/2)])
    z_range = np.int16([peak_coords[2]-np.floor(size_z/2),peak_coords[2]+np.floor(size_z/2)])
    
    mask_bounds = np.int16([np.floor(size_x/2),np.floor(size_y/2), np.floor(size_z/2)])    
    mask_temp = mask
    
    # X-coordinates too small:
    if  peak_coords[0] < mask_bounds[0]:
        x_range = np.int16([0,peak_coords[0]+np.floor(size_x/2)])
        mask_temp = mask_temp[mask_bounds[0]-peak_coords[0]::,:,:]
     
    # X-coordinates too large:
    if peak_coords[0] + mask_bounds[0] >= xpix:
        x_range = np.int16([peak_coords[0]-np.floor(size_x/2),xpix])
        ind_cut = ((x_range[0])+size_x)-xpix
        mask_temp = mask_temp[0:(size_x-ind_cut),:,:]
    
    # Y-coordinates too small:
    if peak_coords[1] < mask_bounds[1]:
        y_range = np.int16([0,peak_coords[1]+np.floor(size_y/2)])
        mask_temp = mask_temp[:,mask_bounds[1]-peak_coords[1]::,:]
        
    # Y-coordinates too large:             
    if peak_coords[1] + mask_bounds[1] >= ypix:
        y_range = np.int16([peak_coords[1]-np.floor(size_y/2),ypix])
        ind_cut = ((y_range[0])+size_y)-ypix
        mask_temp = mask_temp[:,0:(size_y-ind_cut),:]
        
    # Z-coordinates too small:          
    if peak_coords[2] < mask_bounds[2]:
        z_range = np.int16([0,peak_coords[2]+np.floor(size_z/2)])
        mask_temp = mask_temp[:,:,mask_bounds[2]-peak_coords[2]::]
        
    # Z-coordinates too large:
    if peak_coords[2] + mask_bounds[2] >= zsteps:
        z_range = np.int16([peak_coords[2]-np.floor(size_z/2),zsteps])
        ind_cut = ((z_range[0])+size_z)-zsteps
        mask_temp = mask_temp[:,:,0:(size_z-ind_cut)]
        
    ROI_stack = stack[x_range[0]:x_range[1]+1,y_range[0]:y_range[1]+1,z_range[0]:z_range[1]+1]
    mask_size = np.shape(mask_temp)
    
    return ROI_stack, mask_size
    
    
#%% Get Trace Coordinates (X,Y,Z) from a stack
    
def Get_Traces_3D(filepath,stack,threshold_factor,mask_size,max_num_peaks):
        
    # Get file parameters     
    logfile = file_io.read_logfile(filepath)

    ypix = np.uint(str.split(logfile[8],'=')[1])
    xpix = np.uint(str.split(logfile[9],'=')[1])
    zsteps = logfile[18]
    zsteps = str.split(zsteps,"=")[1]
    zsteps = np.uint(str.split(zsteps,",")[0])    
        
    # Make 3D Mask
    size_x, size_y, size_z = mask_size 
    mask_bounds = np.int16([np.floor(size_x/2),np.floor(size_y/2), np.floor(size_z/2)])
    mask = create_3D_mask(size_x, size_y, size_z,10,10)
    
    # Find maximum intensity indeces
    
    threshold = np.median(stack)*threshold_factor
    max_intensity=np.max(stack)
    masked_stack = stack
    peak_coordinates = np.array([], dtype=np.uint16)
    num_peaks = 0
    
    if not max_num_peaks or max_num_peaks is 0:
        max_num_peaks = np.inf
    
    while max_intensity > threshold and num_peaks <= max_num_peaks:
    
        peak_coords = np.unravel_index(np.argmax(masked_stack), np.shape(masked_stack))
        max_intensity = masked_stack[peak_coords]
        peak_coords = np.int16(peak_coords)    
        
        x_range = np.int16([peak_coords[0]-np.floor(size_x/2),peak_coords[0]+np.floor(size_x/2)])
        y_range = np.int16([peak_coords[1]-np.floor(size_y/2),peak_coords[1]+np.floor(size_y/2)])
        z_range = np.int16([peak_coords[2]-np.floor(size_z/2),peak_coords[2]+np.floor(size_z/2)])
        
        mask_temp = mask
    
        # X-coordinates too small:
        if  peak_coords[0] < mask_bounds[0]:
            x_range = np.int16([0,peak_coords[0]+np.floor(size_x/2)])
            mask_temp = mask_temp[mask_bounds[0]-peak_coords[0]::,:,:]
         
        # X-coordinates too large:
        if peak_coords[0] + mask_bounds[0] >= xpix:
            x_range = np.int16([peak_coords[0]-np.floor(size_x/2),xpix])
            ind_cut = ((x_range[0])+size_x)-xpix
            mask_temp = mask_temp[0:(size_x-ind_cut),:,:]
    
        # Y-coordinates too small:
        if peak_coords[1] < mask_bounds[1]:
            y_range = np.int16([0,peak_coords[1]+np.floor(size_y/2)])
            mask_temp = mask_temp[:,mask_bounds[1]-peak_coords[1]::,:]
            
        # Y-coordinates too large:             
        if peak_coords[1] + mask_bounds[1] >= ypix:
            y_range = np.int16([peak_coords[1]-np.floor(size_y/2),ypix])
            ind_cut = ((y_range[0])+size_y)-ypix
            mask_temp = mask_temp[:,0:(size_y-ind_cut),:]
            
        # Z-coordinates too small:          
        if peak_coords[2] < mask_bounds[2]:
            z_range = np.int16([0,peak_coords[2]+np.floor(size_z/2)])
            mask_temp = mask_temp[:,:,mask_bounds[2]-peak_coords[2]::]
            
        # Z-coordinates too large:
        if peak_coords[2] + mask_bounds[2] >= zsteps:
            z_range = np.int16([peak_coords[2]-np.floor(size_z/2),zsteps])
            ind_cut = ((z_range[0])+size_z)-zsteps
            mask_temp = mask_temp[:,:,0:(size_z-ind_cut)]
               
        #z_range=np.int16(z_range), x_range=np.int16(x_range), y_range=np.int16(y_range)
        ROI_stack = stack[x_range[0]:x_range[1]+1,y_range[0]:y_range[1]+1,z_range[0]:z_range[1]+1]
        ROI_stack_masked = ROI_stack*mask_temp
        masked_stack[x_range[0]:x_range[1]+1,y_range[0]:y_range[1]+1,z_range[0]:z_range[1]+1]=ROI_stack_masked
        peak_coordinates = np.append(peak_coordinates,peak_coords)
        
        num_peaks = num_peaks+1
    
    peak_coordinates = np.reshape(peak_coordinates,[3,int(len(peak_coordinates)/3)],1)
    peak_coordinates = np.transpose(peak_coordinates)
    
    return peak_coordinates, num_peaks
    
    
#%% Get the global coordinates of traces from a 3D stack of images 
    
def get_global_coords(filepath, threshold_factor, mask_size,max_num_peaks):
    # Get file parameters 
    logfile = file_io.read_logfile(filepath)
    
    ypix   = np.uint(str.split(logfile[8],'=')[1])
    xpix   = np.uint(str.split(logfile[9],'=')[1])
    zsteps = logfile[18]
    zsteps = str.split(zsteps,"=")[1]
    zsteps = np.uint(str.split(zsteps,",")[0])
    nframes= np.uint(str.split(logfile[7],'=')[1])
    nstacks= int(nframes/zsteps)
    
    # Get Peak coordinates from all stacks 
    stack=np.zeros([xpix,ypix,zsteps], dtype = np.uint16)
    num_traces_pstack=np.empty(nstacks, dtype = np.uint16)
    peak_coordinates = list(range(0,nstacks))
    
    threshold_factor=1.8
    mask_size = [11,11,17]
    for stack_nr in range(0,nstacks):
        # Read stack from file     
        for slice_nr in range(stack_nr*zsteps,stack_nr*zsteps + zsteps):
            stack[:,:,slice_nr-stack_nr*zsteps]=file_io.read_bin(filepath,slice_nr)
            
        # Get Peak coordinates from stack
        [peak_coordinates_stack, num_trace] = Get_Traces_3D(filepath,stack,threshold_factor,mask_size,max_num_peaks)
        peak_coordinates[stack_nr] = peak_coordinates_stack
        num_traces_pstack[stack_nr]=num_trace
    
    return peak_coordinates, num_traces_pstack
    
    
#%% Model function to be used to fit to the data:
    
def gauss_3D(xyz, *p):
    xx,yy,zz = xyz[0],xyz[1],xyz[2]
    A,C,x0,y0,z0,wx,wy,wz = p
    gauss_x = np.square(xx-x0)/(2*np.square(wx))
    gauss_y = np.square(yy-y0)/(2*np.square(wy))
    gauss_z = np.square(zz-z0)/(2*np.square(wz))
    
    return A*np.exp(-1*(gauss_x+gauss_y+gauss_z))+C  
    
#%% Get Local Coordinates by fitting 3D Gaussian to location from Global Coordinates
    
def get_local_coords(filepath,global_coords,mask_size):

    # Get Measurement Parameters 
    logfile = file_io.read_logfile(filepath)
    zsteps = logfile[18]
    zsteps = str.split(zsteps,"=")[1]
    zsteps = np.uint(str.split(zsteps,",")[0])
    
### Get ROI from the stack using Global Coordinates and Fit 3D Gaussian
    
    # Make 3D Mask
    size_x,size_y,size_z=mask_size
    mask = create_3D_mask(size_x,size_y,size_z,10,10)
    
    # Define Bounds
    max_bounds=(65536,65536,size_x,size_y,size_z,size_x-1,size_y-1,size_z-1)
    min_bounds=(0,0,0,0,0,0,0,0)
    bounds = (min_bounds,max_bounds)

### Fit ROI to a 3D Gauss 
      
    # Allocate memory for all fit parameters and errors 
    num_stacks= len(global_coords)    
    fit_params_all, mask_size_all, fit_errors_all = [list(range(0,num_stacks)),list(range(0,num_stacks)),list(range(0,num_stacks))]
    
    for stack_nr in range(0,num_stacks):
        
        # Read one stack from .bin file 
        stack = file_io.get_stack(filepath,stack_nr)
        global_coords_stack = global_coords[stack_nr]
        
        # Allocate memory for fit parameters for 1 stack 
        coords_stack = global_coords[stack_nr]
        fit_params_stack = np.empty([len(coords_stack[:,0]),len(max_bounds)])
        fit_errors_stack = np.empty([len(coords_stack[:,0]),len(max_bounds)])
        mask_size_stack  = np.empty([len(coords_stack[:,0]),3])
           
        for trace_nr in range(0,len(global_coords_stack)):
            
            # Get ROI
            peak_coords = global_coords_stack[trace_nr,:]
            ROI_stack,mask_size = get_ROI_from_stack(filepath,stack,peak_coords,mask)
            
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
            
            # p0 is the initial guess for the fit coefficients
            p0 = [in_A, in_C, in_x0, in_y0, in_z0, in_wx, in_wy, in_wz]
            
            try:
                coeff, var_matrix = curve_fit(gauss_3D, xyz, data, p0=p0,bounds=bounds, absolute_sigma=True)
                perr = np.sqrt(np.diag(var_matrix))
                
            except RuntimeError:
                print('Error - curve_fit failed')
                coeff = np.empty(len(p0))                
                perr  = np.empty(len(p0))
                
            fit_params_stack[trace_nr,:]= coeff
            fit_errors_stack[trace_nr,:]= perr
            mask_size_stack[trace_nr,:] = mask_size

            
        fit_params_all[stack_nr]= fit_params_stack
        fit_errors_all[stack_nr]= fit_errors_stack
        mask_size_all[stack_nr] = mask_size_stack
        
### Transform Global to Local (sub-pixel) Coordinates
    local_coords = list(range(0,num_stacks))
    for stack_nr in range(0,num_stacks):
        
        fit_x,fit_y,fit_z            = fit_params_all[stack_nr][:,2],fit_params_all[stack_nr][:,3],fit_params_all[stack_nr][:,4]
        global_x,global_y,global_z   = global_coords[stack_nr][:,0],global_coords[stack_nr][:,1],global_coords[stack_nr][:,2]
        mask_size_x, mask_size_y, mask_size_z = mask_size_all[stack_nr][:,0],mask_size_all[stack_nr][:,1],mask_size_all[stack_nr][:,2]
        
        centered_x = fit_x-mask_size_x/2
        centered_y = fit_y-mask_size_y/2
        centered_z = fit_z-mask_size_z/2
        
        local_x,local_y,local_z = global_x+centered_x, global_y+centered_y, global_z+centered_z
        local_coords_temp  = np.swapaxes(np.array([local_x,local_y,local_z],order='C'),0,1)
        local_coords[stack_nr] = local_coords_temp
        
    return local_coords, fit_params_all, fit_errors_all