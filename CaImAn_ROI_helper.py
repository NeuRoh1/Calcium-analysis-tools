#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 12 11:51:44 2019

@author: Chang-Eop
"""

import numpy as np
def exclude_rois(signal,excl_rois):
    """ signal: rows are rois, columns are time points
        excl_rois: list of index of rois that will be excluded. Index starts from 0. 
        
        Examples
        >>>F_new = exclude_rois(F, [0,2,3]) 
        -> exclusion of 1st, 3rd, 4th rois from F
        """
    excl_rois = np.array(excl_rois)-1
    result = np.delete(signal, excl_rois, axis = 0)
    return result

def exclude_rois_map(A, excl_rois):
    """ A: rows are flattend pixels, columns are rois
        excl_rois: list of index of rois that will be excluded. Index starts from 0. 
        
        Examples
        >>>A_new = exclude_rois(A, [0,2,3]) 
        -> exclusion of 1st, 3rd, 4th rois from A
        """
    
    excl_rois = np.array(excl_rois)-1
    result = np.delete(A, excl_rois, axis = 1)
    return result    



def margin_detect(input_image, neighbor_type):
    """neighbor_type = 1: left,right, sup, inf 4 pixels
       neighbor_type = 2: 8 neigborhood pixels
    """
    image = input_image.copy()
    remove_i = []
    remove_j = []
    for i in np.arange(1,image.shape[0]-1):
        for j in np.arange(1,image.shape[1]-1):
            if neighbor_type== 1:
                pad = [image[i-1,j], image[i+1,j],image[i,j-1],image[i,j+1]]
                if pad[0]==pad[1]==pad[2]==pad[3]== image[i,j]:
                    remove_i.append(i), remove_j.append(j)
                    

            if neighbor_type== 2:
                pad = [image[i-1,j], image[i+1,j],image[i,j-1],image[i,j+1], image[i-1,j-1], image[i+1,j+1],image[i+1,j-1],image[i-1,j+1]]
                if pad[0]==pad[1]==pad[2]==pad[3]==pad[4]==pad[5]==pad[6]==pad[7]== image[i,j]:
                    remove_i.append(i), remove_j.append(j)
    

    image[remove_i,remove_j] = 0   
                    
    return image     

