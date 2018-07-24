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
filepath = 'C:\\Measurement Data\\180712 - Non_fluorescent_fish_w_lips_and_Lips_on_glass\\data_007.bin'

#%% Get Global Peak Coordinates from the entire range of stacks in the file 

threshold_factor = 2        #Threshold = threshold*median_stack
mask_size = [11,11,23]      #XYZ in pix 
max_num_peaks = 0

global_coords, num_traces = func.get_global_coords(filepath,threshold_factor,mask_size,max_num_peaks)
print(num_traces)

#%% Get Local Peak Coordinates by Fitting 3D Gaussian

mask_size = [11,11,17]

# Local Coords dims=[stack_nr][trace_nr,parameter]
local_coords, fit_params, fit_errors = func.get_local_coords(filepath,global_coords,mask_size)

#%%
#x_loc = np.ndarray.flatten(x_loc)*0.260
#y_loc = np.ndarray.flatten(y_loc)*0.260
#z_loc = np.ndarray.flatten(z_loc)*0.250

num_stacks = len(local_coords)
for stack_nr in range(0,num_stacks):
    
    x_loc = local_coords[stack_nr][:,0]
    y_loc = local_coords[stack_nr][:,1]
    z_loc = local_coords[stack_nr][:,2]

    fig = plt.figure(1) 
    ax = fig.gca(projection='3d')
    ax.scatter(x_loc,y_loc,z_loc,zdir='z')

plt.axis([0,512,0,512])
plt.xlabel('X-Coordinate [pix]')
plt.ylabel('Y-Coorindate [pix]')

plt.figure(2)
plt.title('Fit Errors')
fit_err_x = np.ndarray.flatten(fit_errors[0][:,2])*260
fit_err_y = np.ndarray.flatten(fit_errors[0][:,3])*260
fit_err_z = np.ndarray.flatten(fit_errors[0][:,4])*250

plt.subplot(3,1,1)
plt.title('X-Coordinate')
plt.hist(fit_err_x,50,[0,10])
plt.subplot(3,1,2)
plt.title('Y-Coordinate')
plt.hist(fit_err_y,50,[0,10])
plt.subplot(3,1,3)
plt.title('Z-Coordinate')
plt.hist(fit_err_z,50,[0,10])
plt.tight_layout()

plt.xlabel('Fit Error [nm]')



#%%
A_fit = np.empty(600)
for stack_nr in range(0,600):
    A_fit_temp=fit_params[stack_nr][0,0]
    A_fit[stack_nr]=A_fit_temp
    
A_fit_err = np.empty(600)
for stack_nr in range(0,600):
    A_fit_err_temp=fit_errors[stack_nr][0,0]
    A_fit_err[stack_nr]=A_fit_temp  
    
plt.figure(3)
plt.subplot(1,2,2)
plt.hist(A_fit_err)
plt.subplot(1,2,1)
plt.plot(A_fit)



#%% Nearest Neighbour Method 
# Get trace from stack 1 

stacks_sorted = list(range(0,len(local_coords)-1))
stacks_sorted[0]=local_coords[0]

num_stacks = len(local_coords)

local_coords_temp = list(local_coords)

for stack_nr in range(0,num_stacks-1):
    coord_stack=np.array(stacks_sorted[stack_nr])
    coord_stack_temp = np.array(stacks_sorted[stack_nr])
    
    num_traces = len(local_coords[stack_nr])
    if len(coord_stack) < num_traces:
        num_traces = len(coord_stack)
                
    stack_sorted = np.empty([num_traces,3])

    for trace_nr in range(0,num_traces):
        coord_trace=local_coords[stack_nr][trace_nr]
        r_all = np.sqrt((coord_trace[0]-coord_stack[:,0])**2+(coord_trace[1]-coord_stack[:,1])**2+(coord_trace[2]-coord_stack[:,2])**2)
        
        if len(r_all) is 0:
            continue
        
        else:
            r_min,indices_min = func.get_min(r_all)
            stack_sorted[trace_nr] = coord_stack[indices_min]
            coord_stack_temp[trace_nr] = np.inf
    
    stacks_sorted[stack_nr] = stack_sorted
    
    
#%%
    
plt.figure(4)
for i in range(0,len(stacks_sorted)-1):
    try:
        x,y,z = stacks_sorted[i][0,0],stacks_sorted[i][0,1],stacks_sorted[i][0,2]
        fig = plt.figure(4) 
        ax = fig.gca(projection='3d')
        ax.scatter(x,y,z,zdir='z')
        
    except IndexError:
        print(':-(')
        
            