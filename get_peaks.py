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
from mpl_toolkits.mplot3d import Axes3D

#%% Build Stack from 2D slices 
#Get file path 
filepath = 'C:\\Measurement Data\\180618 - PSF GNRs & DOE Pattern\\data_007.bin'

#%% Get Global Peak Coordinates from the entire range of stacks in the file 
threshold_factor = 2        #Threshold = threshold*median_stack
mask_size = [11,11,23]   #XYZ in pix 
max_num_peaks = 2

global_coords, num_traces = func.get_global_coords(filepath,threshold_factor,mask_size,max_num_peaks)
print(num_traces)

#%% Get Local Peak Coordinates by Fitting 3D Gaussian
mask_size = [11,11,17]
local_coords, fit_params, fit_errors = func.get_local_coords(filepath,global_coords,mask_size)

#%% Get RMS from coordinate position
#x_coordinates 
x_loc = local_coords[0][0]
y_loc = local_coords[1][0]
z_loc = local_coords[2][0]
#%%
x_loc = np.ndarray.flatten(x_loc)*0.260
y_loc = np.ndarray.flatten(y_loc)*0.260
z_loc = np.ndarray.flatten(z_loc)*0.250

fig = plt.figure(1) 
ax = fig.gca(projection='3d')
ax.scatter(x_loc,y_loc,z_loc,zdir='z')
plt.axis([25,225,25,225])

plt.figure(2)
plt.title('Fit Errors')
fit_err_x = np.ndarray.flatten(fit_errors[0,2,:])*260
fit_err_y = np.ndarray.flatten(fit_errors[0,3,:])*260
fit_err_z = np.ndarray.flatten(fit_errors[0,4,:])*250

plt.subplot(3,1,1)
plt.title('X-Coordinate')
plt.hist(fit_err_x,50,[0,10])
plt.subplot(3,1,2)
plt.title('Y-Coordinate')
plt.hist(fit_err_y,50,[0,10])
plt.subplot(3,1,3)
plt.title('Z-Coordinate')
plt.hist(fit_err_z,50,[0,10])

plt.xlabel('Fit Error [nm]')

#%% Plot Intensity Decrease 

fit_param_I = fit_params[0,0,0:380]
plt.figure(3)
plt.plot(fit_param_I)

#%%

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

num_traces=np.empty(nstacks, dtype = np.uint16)
peak_coordinates = np.array([], dtype = np.uint16)
threshold = np.median(stack)*threshold_factor

for stack_nr in range(0,nstacks):
    # Read stack from file     
    for slice_nr in range(stack_nr*zsteps,stack_nr*zsteps + zsteps):
        stack[:,:,slice_nr-stack_nr*zsteps]=file_io.read_bin(filepath,slice_nr)
        
    # Get Peak coordinates from stack
    [peak_coordinates_stack, num_trace] = func.Get_Traces_3D(filepath,stack,threshold,mask_size,max_num_peaks)
    peak_coordinates = np.append(peak_coordinates,peak_coordinates_stack)
    num_traces[stack_nr]=num_trace
    
# Order the found peak coordinates in 3D array 
max_ntraces = int(max(num_traces))

peak_coordinates_global = np.zeros([max_ntraces,3,nstacks], dtype = np.int16)
peak_coordinates = np.reshape(peak_coordinates,[3,int(len(peak_coordinates)/3)],1)
peak_coordinates = np.transpose(peak_coordinates)

for stack_nr in range(0,nstacks):
    peaks_1slice = peak_coordinates[0:num_traces[stack_nr],:]
    peak_coordinates_global[0:len(peaks_1slice),:,stack_nr]=peaks_1slice
    peak_coordinates=peak_coordinates[0:num_traces[stack_nr],:]
        
    