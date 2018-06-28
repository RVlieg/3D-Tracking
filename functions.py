# -*- coding: utf-8 -*-
"""
Created on Thu Jun 28 14:57:53 2018

@author: vlieg
"""

#%% Change Extension of a String

def ChangeExtension(path_file, extension):
    path_old = path_file
    
    extension_index = str.find(path_old, '.')
    
    if str.find(extension,'.') != 0:
        extension = '.' + extension
        
    path_new = path_old[0:extension_index] + extension

    return path_new
    
    


#%% 