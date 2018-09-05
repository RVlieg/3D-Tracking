# -*- coding: utf-8 -*-
"""
Created on Mon Jul  2 14:52:00 2018

@author: vlieg
"""

import numpy as np
import functions as func
import file_io as file_io
import matplotlib.pyplot as plt
import pandas as pd
import func_plot as func_plot


from scipy.optimize import curve_fit
from mpl_toolkits.mplot3d import Axes3D

#%% Build Stack from 2D slices 
#Get file path 
filepath = file_io.get_filepath()

#%% Get Global Peak Coordinates from the entire range of stacks in the file 

threshold_factor = 2        #Threshold = threshold*median_stack
mask_size = [11,11,23]      #XYZ in pix 
max_num_peaks = 0

global_coords, num_traces = func.get_global_coords(filepath,threshold_factor,mask_size,max_num_peaks)
print(num_traces)


#%% Write Global coordinates (= list) to .xlsx file    
file_io.write_xlsx_list(filepath,global_coords,['x [pix]','y [pix]','z [pix]'],[0,1,2])

#%% Get Local Peak Coordinates by Fitting 3D Gaussian
mask_size = [11,11,17]

# Local Coords dims=[stack_nr][trace_nr,parameter]
local_coords, fit_params, fit_errors = func.get_local_coords(filepath,global_coords,mask_size)

#%% Write the output of the 3D Gaussian fit to .xlsx file
file_io.write_xlsx_list(filepath,local_coords,['x_local [pix]','y_local [pix]','z_local [pix]'],[3,4,5])
file_io.write_xlsx_list(filepath,fit_params,['A [a.u]','C [a.u.]','x0 [pix]','y0 [pix]','z0 [pix]','wx [pix]','wy [pix]','wz [pix]'],[6,7,8,9,10,11,12,13])
file_io.write_xlsx_list(filepath,fit_errors,['A_err [a.u]','C_err [a.u.]','x_err [pix]','y0_err [pix]','z0_err [pix]','wx_err [pix]','wy_err [pix]','wz_err [pix]'],[14,15,16,17,18,19,20,21])


#%% Nearest Neighbour Method 

# Sort Stacks to get correct trace numbers per trace based on nearest distance (r) 


stacks_sorted = list(range(0,len(local_coords)))
stacks_sorted[0]=local_coords[0]

num_stacks = len(local_coords)

local_coords_temp = list(local_coords)
trace_nrs = list(range(0,len(local_coords)))

trace_nrs[0]=np.arange(0,np.shape(local_coords[0])[0])
r_all_stacks = list(range(0,len(local_coords)))
r_all_stacks[0] = np.empty([len(local_coords[0]),len(local_coords[0])])

for stack_nr in range(0,len(global_coords)-1):
#for stack_nr in range(0,3):
    #Allocate memory for coordinates of second stack 
    stack2=np.array(local_coords[stack_nr+1])
    
    
    #Make temporary copy so values can be edited during iteration
    coord_stack_temp = np.array(local_coords[stack_nr+1])
    
    num_traces = len(local_coords[stack_nr])
    
    #Allocate memory for all the sorted trace numbers and distances (r)
    trace_nrs[stack_nr+1] = np.empty(len(stack2),dtype=np.int)
    r_all_stacks[stack_nr+1] = np.empty([len(stack2),len(stack2)])

    
    if len(stack2) < len(local_coords[stack_nr]):
        num_traces = len(stack2)

    stack_sorted = np.empty([num_traces,3])
    
    
    for trace_nr in range(0,num_traces):
        
        #coord_trace = Trace from first stack which is compared with second stack 
        stack1_trace = local_coords[stack_nr][trace_nr,:]
        
        r_all = np.sqrt((stack1_trace[0]-coord_stack_temp[:,0])**2+(stack1_trace[1]-coord_stack_temp[:,1])**2+(stack1_trace[2]-coord_stack_temp[:,2])**2)
        r_all_stacks[stack_nr+1][:,trace_nr] = r_all
        
        if len(r_all) is 0:
            stack_sorted = np.empty([num_traces,3])
        
        else:
            #Get minimum r-value and corresponding index of second stack 
            r_min, indices_min = func.get_min(r_all)
            
            #Sort matched trace of stack2 to same index as stack1 
            trace_couple = stack2[indices_min]
            stack_sorted[trace_nr] = trace_couple
            coord_stack_temp[indices_min] = np.inf
            
            #Link the found index_rmin to the corresponding trace in stack2
            trace_nr_rmin = trace_nrs[stack_nr][trace_nr]
            trace_nrs[stack_nr+1][indices_min] = trace_nr_rmin
            
        #When length of stack2 is larger than stack1, append numbers starting from max(index_number_stack1) after all iterations
        if trace_nr == num_traces-1 and len(stack2) > len(local_coords[stack_nr]):
            for i in range(0,len(trace_nrs[stack_nr+1])):
                if trace_nrs[i,:] not []:
                    
        
            #trace_nrs[stack_nr+1][trace_nr+1::] = [max(trace_nrs[stack_nr])+1:len(stack_nr+1)]       
            
        #When length of stack2 is  than stack1, append numbers starting from max(index_number_stack1) after all iterations
        #elif trace_nr == num_traces-1 and len(stack2) < len(local_coords[stack_nr]):
         #   trace_nrs[stack_nr+1][trace_nr+1::]=np.arange(max(trace_nrs[stack_nr])+1,len(trace_nrs[stack_nr+1]))
            
    # Attach all left-over coordinates to the sorted stack when stack2 is larger than stack1  
    if len(stack2) > num_traces:
        stack_sorted = np.append(stack_sorted,stack2[num_traces+1::],0) 
        
    stacks_sorted[stack_nr+1] = stack_sorted



#%% Plot Data

func_plot.plot_hist_errs(fit_errors,1,filepath)



