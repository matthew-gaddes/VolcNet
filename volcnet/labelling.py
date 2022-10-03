#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 10 11:47:46 2022

@author: matthew


"""

import pdb

#%%


def label_volcnet_files(volcnet_files, def_min = 0.05):
    """ Given a list of volcnet files, create all possible labels for the interferograms that can be made between all possible acquisitions.  
    Note that this doesn't actually make the interferograms, so the outputs are small (i.e. not GBs)
    
    Inputs:
        volcnet_files | list of strings | list of volcnet files to label.  
        def_min | float | magnitude of deformation must be larger than this to be classed as deformation.  Units: metres.  
        
    Returns:
        labels_dyke | many x 4 | Interferograms that are labelled as dyke.  First coumn is file number, second is acquisition 1, third is acquisition 2, fourth is deformation magnitude
        labels_sill | many x 4 |
        labels_atmo | many x 4 |

    History:
        2022_10_03 | MEG | Written
    """
    import numpy as np
    import pickle
        
    labels_dyke = np.nan * np.zeros((1, 4))                                         # initialise, 1 x 4.  First coumn is file number, second is acquisition 1, third is acquisition 2, fourth is deformation magnitude
    labels_sill = np.nan * np.zeros((1, 4))                                          # as above, for next label type
    labels_atmo = np.nan * np.zeros((1, 4))                                          # as above, for next label type
    
    for file_n, volcnet_file in enumerate(volcnet_files):                               # loop through files.  
        
        print(f"Opening file: {volcnet_file.split('/')[-1]}")
        # 1: Open the file
        with open(volcnet_file, 'rb') as f:
            displacement_r3 = pickle.load(f)
            tbaseline_info = pickle.load(f)
            persistent_defs = pickle.load(f)
            transient_defs = pickle.load(f)
        print(f"The interferograms are of size: {displacement_r3['mask'].shape}")
        #n_acq, ny, nx = displacement_r3['cumulative'].shape
        
       
        for acq_n1, acq_1 in enumerate(tbaseline_info['acq_dates']):                                                # 
            for acq_n2, acq_2 in enumerate(tbaseline_info['acq_dates']):
                if acq_1 == acq_2:                                                                                                          # will just be zeros so ignore.   
                    pass                                                                                                                      # just leave as nans          
                else:
                    #ifg = displacement_r3['cumulative'][acq_n2,] - displacement_r3['cumulative'][acq_n1,]                                    # make the ifg between the two acquisitions.  
                    def_predicted, sources, def_location = label_volcnet_ifg(f"{acq_1}_{acq_2}", persistent_defs, transient_defs)              # label the ifg, def_location is still in terms of lon and lat
                    if (np.abs(def_predicted) < def_min):                                                                                    # if the deformation is less than the threshold selected       
                        # def_magnitudes[acq_n1, acq_n2] = 0.                                                                                  # set deformation to 0.           
                        # labels[acq_n1, acq_n2] = 2                                                                                           # and label as no def
                        labels_atmo[-1] = np.array([file_n, acq_n1, acq_n2, def_predicted])
                        labels_atmo = np.vstack((labels_atmo, np.nan * np.zeros((1, 4))))
                    else:                                                                                                                    # else deformation is big enough
         #               def_magnitudes[acq_n1, acq_n2] = def_predicted                                                                       # record deformation
                        if sources[0] == 'dyke':
                            # labels[acq_n1, acq_n2] = 0                                                                                      # label as dyke                    
                            labels_dyke[-1] = np.array([file_n, acq_n1, acq_n2, def_predicted])
                            labels_dyke = np.vstack((labels_dyke, np.nan * np.zeros((1, 4))))
                        elif sources[0] == 'sill':
    #                        labels[acq_n1, acq_n2] = 1                                                                                      # label as sill
                            labels_sill[-1] = np.array([file_n, acq_n1, acq_n2, def_predicted])
                            labels_sill = np.vstack((labels_sill, np.nan * np.zeros((1, 4))))
        #data_indexs[volcnet_file.split('/')[-1]] = (def_magnitudes, labels)
    labels_dyke = labels_dyke[:-1,:]
    labels_sill = labels_sill[:-1,:]
    labels_atmo = labels_atmo[:-1,:]
    return labels_dyke, labels_sill, labels_atmo




#%%

def label_volcnet_ifg(ifg_name, persistent_defs, transient_defs, noise_threshold = 0.02):
    """ Given VolcNet labels (persistent and transient defs) for a time series, and the acquisiiton dates for an inteferogram,
    Create a label for it, calculate the amount of deformaiton expected, and the bounding box.  
    
    Inputs:
        ifg_name | string | in form yyyymmdd_yyyymmdd
        persistent_defs | list of dicts | Info on each type of persistent deformation in the time series.  
        transient_defs | list of dicts | Info on each type of transient deformation in the time series.  
        
    Returns:
        def_ifg | float | Amount of deforamtion predicted to be in that interferogram (m)
        sources | list of string | source label for that deformation.  If multiple sources (e.g. a long ifg has both a sill and dyke), list has multiple entries.  
        def_location | list of tuple (lon, lat) | Closed polygon around deformation.  
        
    History:
        2022_05_04 | MEG | Written.  
    
    """
    
    import numpy as np
    from datetime import datetime
    
    from collections import namedtuple
    from shapely.geometry import Polygon
    from shapely.ops import unary_union
        
    def polygon_to_list_tuples(polygon):
        """Given a shapely polygon, turn it back into a simple list of tuples.  
        """
        polygon_list = []
        
        x, y = polygon.exterior.coords.xy                                 # get the coords of that
        x = np.array(x)                                                                     # make into a numpy array, rank 1
        y = np.array(y)                                                                     # make into a numpy array, rank 1
        
        for lon, lat in zip(x, y):
            polygon_list.append((lon, lat))
        return polygon_list
    

    
    Range = namedtuple('Range', ['start', 'end'])
    
    acq_start_dt = datetime.strptime(ifg_name[:8], '%Y%m%d')
    acq_stop_dt = datetime.strptime(ifg_name[9:], '%Y%m%d')
    if acq_start_dt > acq_stop_dt:
        acq_start_dt = datetime.strptime(ifg_name[9:], '%Y%m%d')                # flip so now a standard ifg  1st date -> 2nd date (where 1st comes before)
        acq_stop_dt = datetime.strptime(ifg_name[:8], '%Y%m%d')                 # continue flip
        backward_ifg = True
    else:
        backward_ifg = False
        
    tbaseline = (acq_stop_dt - acq_start_dt).days
    if backward_ifg:                                                                # if a backward ifg
        tbaseline *= (-1)                                                           # t baseline is negative  
    range_ifg = Range(start = acq_start_dt, end = acq_stop_dt)
    
    def_ifg = 0.                                                                                      # initiate
    sources = []                                                                                            # initiate.  
    location_polygon = Polygon([])
    
    # 1 Add any persistent deformation
    for persistent_def in persistent_defs:
        d_start = datetime.strptime(str(persistent_def['def_episode_start']), '%Y%m%d')
        d_stop = datetime.strptime(str(persistent_def['def_episode_stop']), '%Y%m%d')
        range_def = Range(start = d_start, end = d_stop)
        latest_start = max(range_ifg.start, range_def.start)
        earliest_end = min(range_ifg.end, range_def.end)
        delta = (earliest_end - latest_start).days
        overlap = max(0, delta)         
        # print(f"Overlapping days: {overlap}")                                                                     # can be useful to debug/test
        if overlap > 0:
            def_episode = ((overlap / 365.25) * persistent_def['def_rate'])                                      # convert days to years, then multiply by rate in m/year, to get deforamtion in the ifg due to this persistent episode.  
            if backward_ifg:                                                                                     # but if it's a backward ifg...
                def_episode *= (-1)                                                                             # signal will be in opposite sense
            def_ifg += def_episode
            if persistent_def['source'] not in sources:
                sources.append(persistent_def['source'])
            location_polygon_current = Polygon(persistent_def['def_polygon'])
            location_polygon = unary_union([location_polygon, location_polygon_current])
                        
    #print(f"{def_ifg} m")
    
    # 2: Add any transient deformation
    for transient_def in transient_defs:
        d_start = datetime.strptime(str(transient_def['def_episode_start']), '%Y%m%d')
        d_stop = datetime.strptime(str(transient_def['def_episode_stop']), '%Y%m%d')
        range_def = Range(start = d_start, end = d_stop)
        latest_start = max(range_ifg.start, range_def.start)
        earliest_end = min(range_ifg.end, range_def.end)
        delta = (earliest_end - latest_start).days
        overlap = max(0, delta)
        #print(f"Overlapping days: {overlap}")                                                                      # can be useful to debug/test
        if overlap > 0:
            def_episode = transient_def['def_magnitude']                                                        # 
            if backward_ifg:                                                                                    # but if it's a backward ifg...
                def_episode *= (-1)                                                                             # signal will be in opposite sense
            def_ifg += def_episode
            if transient_def['source'] not in sources:
                sources.append(transient_def['source'])
            location_polygon_current = Polygon(transient_def['def_polygon'])
            location_polygon = unary_union([location_polygon, location_polygon_current])
            
    def_location = polygon_to_list_tuples(location_polygon)
    
    
    return def_ifg, sources, def_location

#%%