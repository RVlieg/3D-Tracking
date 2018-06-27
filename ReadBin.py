# -*- coding: utf-8 -*-
"""
Created on Thu Mar 15 11:15:31 2018

@author: Redmar Vlieg
"""

# Import all the necessary Python Packages 
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import tkinter as tk
from tkinter import filedialog
import array
import struct

#%% Get the Filename
root = tk.Tk()
root.withdraw()

filename = filedialog.askopenfilename()
#filename = filename.replace('/','\\')


#%% Open the File 
f = open(filename, 'rb')
data = np.fromfile(f, dtype = np.uint16)
nz = np.int(len(data)/(400*400))
data = np.reshape(data,[400,400,nz],'F')
f.close()
plt.imshow(data[:,:,5])

data2 = np.array([])
with open(filename,'rb') as f:
    byte_code=f.read(1)
    while byte_code != "":
        byte=f.read(1)
        data2=np.append(data2,byte)
        
data3 = np.asarray(data2, dtype=np.uint16)

