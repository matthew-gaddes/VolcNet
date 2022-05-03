#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  3 14:34:41 2022

@author: matthew
"""

import os
import glob
from pathlib import Path
from copy import deepcopy
import pdb

def read_volcnet_label(label_file):
    """Given a .txt file of labels (i.e. deformation or no deformation etc), read it into a dictionary.  
    Inputs:
        label_file | pathlib Path | .txt file to open
    Returns:
        volcnet_labels
       
    History:
        2022_04_20 | MEG | Modify from similar function in LiCSAlert
    """
    import configparser    
   
    volcnet_labels = {}
   
    config = configparser.ConfigParser()                                                       # read the config file
    config.read_file(open(label_file))
    
    persistent_deformation_events = []
    transient_deformation_events = []
    
    for deformation_event in config.sections():
        event = {}
        
        # first open the parameters that are shared for both persistent and transient deformation
        event['def_lon_west'] = float(config.get(deformation_event, 'def_lon_west'))       # 
        event['def_lon_east'] = float(config.get(deformation_event, 'def_lon_east'))            # 
        event['def_lat_south'] = float(config.get(deformation_event, 'def_lat_south'))       # 
        event['def_lat_north'] = float(config.get(deformation_event, 'def_lat_north'))       # 
        event['def_episode_start'] = int(config.get(deformation_event, 'def_episode_start'))       # 
        event['def_episode_stop'] = int(config.get(deformation_event, 'def_episode_stop'))       # 
        event['source'] = config.get(deformation_event, 'source')       # 
        
        if deformation_event[:10] == 'persistent':                                                          # then open only persistent event paramters
            event['def_rate'] = float(config.get(deformation_event, 'def_rate'))       # 
            persistent_deformation_events.append(event)
            
        else:                                                                                               # or transient event parameters.  
            event['def_magnitude'] = float(config.get(deformation_event, 'def_magnitude'))       # 
            transient_deformation_events.append(event)
        
    return persistent_deformation_events, transient_deformation_events


def bounding_to_polygon(dict_bounding):
    """
    """
    from copy import deepcopy
    
    dict_polygon = deepcopy(dict_bounding)
    del dict_polygon['def_lon_west'], dict_polygon['def_lon_east'], dict_polygon['def_lat_south'], dict_polygon['def_lat_north']
    
    west = dict_bounding['def_lon_west']
    east = dict_bounding['def_lon_east']
    south = dict_bounding['def_lat_south']
    north = dict_bounding['def_lat_north']
    polygon = f"[({west:.2f}, {north:.2f}), ({east:.2f}, {north:.2f}), ({east:.2f}, {south:.2f}), ({west:.2f}, {south:.2f}), ({west:.2f}, {north:.2f})]"                # clockwise, closed, preserve 2 dp
    dict_polygon['def_polygon'] = polygon
    

    return dict_polygon

#%%

annotation_dirs = [Path("/home/matthew/university_work/13_volcnet/raw_annotation_data/fabien_agung"),
                   Path("/home/matthew/university_work/13_volcnet/raw_annotation_data/licsbas"),
                   Path("/home/matthew/university_work/13_volcnet/raw_annotation_data/marco_galapagos")]


outdir = Path("/home/matthew/university_work/13_volcnet/raw_annotation_data2")

#%%


for annotation_dir in annotation_dirs:
    
    label_files = sorted(glob.glob(str(annotation_dir/ '*.txt')), key = os.path.getmtime)            # get the paths to the mat files from fabien
    
    for label_file in label_files:
        
        persistent_defs, transient_defs = read_volcnet_label(label_file)                            # read the file
        
        persistent_defs_new = []
        for persistent_def in persistent_defs:
            persistent_defs_new.append(bounding_to_polygon(persistent_def))
            
        transient_defs_new = []
        for transient_def in transient_defs:
            transient_defs_new.append(bounding_to_polygon(transient_def))
            
        
        f = open(outdir / label_file.split('/')[-2] / label_file.split('/')[-1], 'w')
        
        for n, persistent_def_new in enumerate(persistent_defs_new):
            f.write(f"[persistent_deformation_{n:02d}]\n")
            for key, value in persistent_def_new.items():
                f.write(f"{key} = {value}\n")
            f.write("\n")
        
        
        for n, transient_def_new in enumerate(transient_defs_new):
            f.write(f"[transient_deformation_{n:02d}]\n")
            for key, value in transient_def_new.items():
                f.write(f"{key} = {value}\n")
            f.write("\n")
            

        f.close()
            
        
        
        
        
        
        