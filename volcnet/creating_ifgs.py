#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  5 15:00:52 2022

@author: matthew
"""


def create_volcnet_ifgs(labels_all, volcnet_files, outdir, n_data_per_file = 100, ny = 224, nx = 224, volcnet_def_min = 0.05):
    """ Given the file numbers, acquisition 1 and acqustion 2 numbers in labels all, and the volcnet files, make those interferograms.  
    Inputs:
        labels_all | many x 4 |  First coumn is file number, second is acquisition 1, third is acquisition 2, fourth is deformation magnitude
        volcnet_files | list of strings | list of volcnet files to label.  
        def_min | float | magnitude of deformation must be larger than this to be classed as deformation.  Units: metres.  
        outdir | pathlib Path |  out directory.  
        n_data_per_file | int | number of data per .pkl file
        ny | int | output size, in pixels
        nx | int | as above. 
        volcnet_def_min | float | deformation must be above this (in metres) be classed as deformation.  
        
    Returns:
        .pkl files with X, Y_class, Y_loc
        
    History:
        2022_10_05 | MEG | written.  
        
    """
    import numpy as np
    import numpy.ma as ma
    import pickle
    

    from volcnet.labelling import label_volcnet_ifg
    from volcnet.aux import ll_2_pixel
    
    from deep_learning_tools.data_handling import rescale_timeseries, random_cropping
    
    def initialise_arrays(n_data, ny, nx,n_channels):
        """
        """
        X = ma.zeros((n_data, ny, nx, n_channels))               # initialise, rank 4 ready for Tensorflow, last dimension is being used for different crops.  
        Y_class = np.zeros((n_data, 3))                                                                             # initialise, doesn't need another dim as label is the same regardless of the crop.  
        Y_loc = np.zeros((n_data, 4))                                                                            # initialise
        return X, Y_class, Y_loc
    
    file_n = 0
    data_n = 0
    X, Y_class, Y_loc = initialise_arrays(n_data_per_file, ny, nx, 1)                # initiliase for next file.      
    
    for ifg_n, ifg_info in enumerate(labels_all):
        
        # 0: open the volcnet file and check not already open)
        if ifg_n == 0:
            current_file = int(ifg_info[0])
            open_volcnet_file = True
    
        else:
            if current_file == int(ifg_info[0]):
                open_volcnet_file = False
            else:
                current_file = int(ifg_info[0])
                open_volcnet_file = True
        
        if open_volcnet_file:
            print(f"Opening file: {volcnet_files[current_file].split('/')[-1]}")
            with open(volcnet_files[current_file], 'rb') as f:
                displacement_r3 = pickle.load(f)
                tbaseline_info = pickle.load(f)
                persistent_defs = pickle.load(f)
                transient_defs = pickle.load(f)
                
            n_acq, ny_original, nx_original = displacement_r3['cumulative'].shape
            n_ifg = (n_acq*n_acq) - n_acq
            print(f"The interferograms are of size: {displacement_r3['mask'].shape}")
           
            if (nx_original < nx) or (ny_original < ny):
                if ny_original/nx_original < 1:                                                                                                       # this is less than 1 if the image is wider than tall.  
                    rescale_factor = (1.4 * ny) /ny_original                                                   #rescale to ensure y is large enough to be cropped down to 224    
                else:
                    rescale_factor = (1.4 * ny) / nx_original
                    
                displacement_r3 = rescale_timeseries(displacement_r3, rescale_factor)
                print(f"The interferograms have been interpolated to size: {displacement_r3['mask'].shape}")
            
         # 1:
        acq_n1 = int(ifg_info[1])                                                                                                   # stored as a float, so convert
        acq_n2 = int(ifg_info[2])
        acq_1 = tbaseline_info['acq_dates'][acq_n1]                                                                                 # get the acquisition date in form YYYYMMDD instead of a number
        acq_2 = tbaseline_info['acq_dates'][acq_n2]
        
        ifg = displacement_r3['cumulative'][acq_n2,] - displacement_r3['cumulative'][acq_n1,]                                     # make the ifg between the two acquisitions.  
        def_predicted, sources, def_location = label_volcnet_ifg(f"{acq_1}_{acq_2}", persistent_defs, transient_defs)             # label the ifg, def_location is still in terms of lon and lat
        def_loc_pixels = ll_2_pixel(def_location, displacement_r3['lons'], displacement_r3['lats'])                               # convert the location label from lon lat to pixels (x then y)
    
        if (np.abs(def_predicted) < volcnet_def_min):                                                                             # if the deformation is less than the threshold selected       
            ifg_cropped_r3 = random_cropping(ifg, ny, None)    
        else:                                                                                                                       # else deformation is big enough
            ifg_cropped_r3, Y_loc_cropped = random_cropping(ifg, ny, def_loc_pixels)    
        
    
        for crop_n in range(9):
            if data_n < (n_data_per_file - 1):                                                                                   # check that aren't overfilling X
                X[data_n,] = ifg_cropped_r3[:,:,crop_n : crop_n + 1]
                if (np.abs(def_predicted) < volcnet_def_min):                                                                              # if the deformation is less than the threshold selected       
                    Y_class[data_n,] = np.array([0,0,1])                                                                                    # this is the one hot encoding for atmo
                    Y_loc[data_n,] = np.array([0,0,0,0])
                else:                                                                                                                       # else deformation is big enough
                    if sources[0] == 'dyke':
                        Y_class[data_n,] = np.array([1,0,0])                                                                                # this is the one hot encoding for dyke
                    elif sources[0] == 'sill':
                        Y_class[data_n,] = np.array([0,1,0])                                                                                # this is the one hot encoding for sill
                    Y_loc[data_n] = Y_loc_cropped[crop_n, ]
                data_n += 1
            else:
                pass
                
        # When generated the required number per file, save the file.  
        if data_n == (n_data_per_file - 1):                         
            print(f"    Saving file {file_n}")
            with open(outdir / f"data_file_unshuffled_{file_n:05d}.pkl", 'wb') as f:                     # save the output as a pickle
                pickle.dump(X, f)
                pickle.dump(Y_class, f)
                pickle.dump(Y_loc, f)
            file_n += 1                                                                                                                                         # advance to next file
            data_n = 0                                                                                                                                          # initiate for next file        
            X, Y_class, Y_loc = initialise_arrays(n_data_per_file, ny, nx, 1)            # initiliase for next file.  
        
        


