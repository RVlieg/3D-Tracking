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
filepath = 'C:\\Users\\Redmar\\Desktop\\data_003.bin'


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

# Read one stack from .bin file 
stack = file_io.get_stack(filepath,0)

#%% Get ROI from the stack using Global Coordinates and Fit 3D Gaussian

# Make 3D Mask
size_x,size_y,size_z=[11,11,7]
mask = func.create_3D_mask(size_x,size_y,size_z,10,10)

# Get ROI
peak_coords = global_coords[15,:,0]
ROI_stack,mask_size = func.get_ROI_from_stack(filepath,stack,peak_coords,mask)

plt.plot(np.ndarray.flatten(ROI_stack))

##### Fit Data 


# Define some test data which is close to Gaussian
data = np.ndarray.flatten(ROI_stack)
size_x,size_y,size_z = mask_size

x = np.arange(0,size_x)
y = np.arange(0,size_y)
z = np.arange(0,size_z)

# Create xyz-coordinates 
xx, yy, zz = np.meshgrid(x, y, z, sparse=False)
xyz = [np.ndarray.flatten(xx),np.ndarray.flatten(yy), np.ndarray.flatten(zz)]



# Define model function to be used to fit to the data above:
def gauss_3D(xyz, *p):
    xx,yy,zz = xyz[0],xyz[1],xyz[2]
    A,x0,y0,z0,wx,wy,wz = p
    gauss_x = np.square(xx-x0)/(2*np.square(wx))
    gauss_y = np.square(yy-y0)/(2*np.square(wy))
    gauss_z = np.square(zz-z0)/(2*np.square(wz))
    return A*np.exp(-1*(gauss_x+gauss_y+gauss_z))

# Initial Fit Parameters 
in_A  = float(max(data))
in_x0 = max(x)/2
in_y0 = max(y)/2
in_z0 = max(z)/2
in_wx = 1.5
in_wy = 1.5
in_wz = 3

# p0 is the initial guess for the fitting coefficients (A, mu and sigma above)
p0 = [in_A, in_x0, in_y0, in_z0, in_wx, in_wy, in_wz]

coeff, var_matrix = curve_fit(gauss_3D, xyz, data, p0=p0)
coeff = list(coeff)

# Plot Found Coefficients 
fitted_curve = gauss_3D(xyz,*coeff)

plt.plot(fitted_curve)
plt.plot(data)
