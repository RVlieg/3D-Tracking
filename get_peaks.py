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

# Read one stack from .bin file 
stack = file_io.get_stack(filepath,0)

#%% Get ROI from the stack using Global Coordinates and Fit 3D Gaussian

# Make 3D Mask
size_x,size_y,size_z=[11,11,23]
mask = func.create_3D_mask(size_x,size_y,size_z,10,10)

# Get ROI
peak_coords = global_coords[0,:,0]
ROI_stack = func.get_ROI_from_stack(filepath,stack,peak_coords,mask)

plt.plot(np.ndarray.flatten(ROI_stack))

#%% Fit Data 


# Define some test data which is close to Gaussian
data = np.ndarray.flatten(ROI_stack)

x = np.arange(0,size_x)
y = np.arange(0,size_y)
z = np.arange(0,size_z)
mid_indx = np.uint16(np.floor(len(x)/2))
mid_indy = np.uint16(np.floor(len(y)/2))
mid_indz = np.uint16(np.floor(len(z)/2))

# Create xyz-coordinates 
xx, yy, zz = np.meshgrid(x, y, z, sparse=False)

# Define model function to be used to fit to the data above:
def gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

def gauss3D(x, *p):
    
    gauss_x = np.square(xx-x[mid_indx])/(2*np.square(S))
    gauss_y = np.square(yy-y[mid_indy])/(2*np.square(S))
    gauss_z = np.square(zz-z[mid_indz])/(2*np.square(S))
    gauss_xyz = A*np.exp(-1*(gauss_x+gauss_y+gauss_z))

# p0 is the initial guess for the fitting coefficients (A, mu and sigma above)
p0 = [1., 0., 1.]

coeff, var_matrix = curve_fit(gauss, bin_centres, hist, p0=p0)

# Get the fitted curve
hist_fit = gauss(bin_centres, *coeff)

plt.plot(bin_centres, hist, label='Test data')
plt.plot(bin_centres, hist_fit, label='Fitted data')

# Finally, lets get the fitting parameters, i.e. the mean and standard deviation:
print('Fitted mean = ', coeff[1])
print('Fitted standard deviation = ', coeff[2])

plt.show()
