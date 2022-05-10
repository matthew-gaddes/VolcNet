#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 10 16:32:36 2022

@author: matthew
"""

#%%

def ll_2_pixel(lls, lons, lats):
    """
    Given an image that uses matrix notation (i.e. 0,0 is top left), get the pixel number of a lon lat coordinate pair.  
    
    Inputs:
        lls | list of tuples | lon lat pairs for each coordinate
        lons | rank 2 array | lons for each pixel
        lats | rank 2 array | lats for each pixel
        
    Returns:
        xys_array | ? x 2 | x and y pixel number (as floats) for each coordinate paire
        
    History:
        2022_05_10 | MEG | Written.  
    """
    
    import numpy as np
    
    xys = np.zeros((len(lls), 2))
    
    ll_upper_left = np.array([lons[0,0], lats[0,0]])                                  # x then y, lon then lat
    
    pixel_size = {}
    pixel_size['x'] = np.abs(lons[0,0] - lons[0,1])                                 # how many degrees between each pixel, in x
    pixel_size['y'] = np.abs(lats[0,0] - lats[1,0])
    
    for point_n, ll in enumerate(lls):
        xys[point_n, 0] = (ll[0] - ll_upper_left[0])   / pixel_size['x']
        xys[point_n, 1] = ((-1) * (ll[1] - ll_upper_left[1]))   / pixel_size['y']         # y needs flipping as down is positive in matrix notation
        
    return xys


#%%