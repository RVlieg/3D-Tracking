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
threshold = 1.5        #Threshold = threshold*median_stack
mask_size = [11,11,23]   #XYZ in pix 

global_coords, num_traces = func.get_global_coords(filepath,threshold,mask_size)
print(num_traces)

#%% Get Local Peak Coordinates by Fitting 3D Gaussian
mask_size = [11,11,17]
local_coords, fit_params, fit_errors = func.get_local_coords(filepath,global_coords,mask_size)

#%% Get RMS from coordinate position
#x_coordinates 
x_loc = local_coords[0]
y_loc = local_coords[1]
z_loc = local_coords[2]

x_loc = np.ndarray.flatten(x_loc)*260
y_loc = np.ndarray.flatten(y_loc)*260
z_loc = np.ndarray.flatten(z_loc)*250

fig = plt.figure(1) 
ax = fig.gca(projection='3d')
ax.scatter(x_loc,y_loc,z_loc,zdir='z')
plt.axis([25,225,25,225])

plt.figure(2)
plt.title('Fit Errors')
fit_err_x = np.ndarray.flatten(fit_errors[:,2,:])*260
fit_err_y = np.ndarray.flatten(fit_errors[:,3,:])*260
fit_err_z = np.ndarray.flatten(fit_errors[:,4,:])*250

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

