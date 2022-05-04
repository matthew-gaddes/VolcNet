#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 15:10:50 2022

@author: matthew
"""

import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import sys
from pathlib import Path
from copy import deepcopy
import pdb
import pickle
import glob
import os
from datetime import datetime



import volcnet
from volcnet.plotting import volcnet_ts_visualiser


# sys.path.append("/home/matthew/university_work/23_insar_tools")                  # 
# import insar_tools
# from insar_tools.plotting import xticks_every_nmonths

# MEG debug imports
sys.path.append("/home/matthew/university_work/python_stuff/python_scripts")
from small_plot_functions import matrix_show, pngs_to_gif, col_to_ma, quick_linegraph
import matplotlib.pyplot as plt
#plt.switch_backend('Qt5Agg')                                                               # only needed if importing this mid debug 





#%% Settings

volcnet_dir = Path("/home/matthew/university_work/13_volcnet")

figsize = 10

#%% Get a list of Volcnet files




volcnet_files = sorted(glob.glob(str(volcnet_dir / '*.pkl')))            # get the paths to the mat files from fabien

# for volcnet_file in volcnet_files[1:2]:
#     print("TESTING - only using Campi Flegrei volcnet file.  ")
# for volcnet_file in volcnet_files[11:12]:
#     print("TESTING - only using Wolf volcnet file.  ")
for volcnet_file in volcnet_files[10:11]:
    print("TESTING - only using Sierra Negra 128 volcnet file.  ")
    
    # 1: Open the file
    with open(volcnet_file, 'rb') as f:
        displacement_r3 = pickle.load(f)
        tbaseline_info = pickle.load(f)
        persistent_defs = pickle.load(f)
        transient_defs = pickle.load(f)

    #print(volcnet_file)
    volcnet_ts_visualiser(displacement_r3, tbaseline_info, persistent_defs, transient_defs, acq_spacing = 5, ifg_resolution = 20, figsize_height = 10)
    
    
        
#%% create labels

    
# make an ifg.  
# check if transient

# check if persistent def exists.  
# estimate rate from persisent.  
# compare to threshold
# assign label.  

def volcnet_labeller(ifg_name, persistent_defs, transient_defs, noise_threshold = 0.02):
    """ Given VolcNet lables (persistent and transient defs) for a time series, and the acquisiiton dates for an inteferogram,
    Create a label for it.  
    
    Inputs:
        
    Returns:
    
    """
    from collections import namedtuple
    
    Range = namedtuple('Range', ['start', 'end'])
    
    acq_start_dt = datetime.strptime(ifg_name[:8], '%Y%m%d')
    acq_stop_dt = datetime.strptime(ifg_name[9:], '%Y%m%d')
    tbaseline = (acq_stop_dt - acq_start_dt).days
    r1 = Range(start = acq_start_dt, end = acq_stop_dt)
    
    def_predicted = 0.                                                                                      # initiate
    sources = []                                                                                            # initiate.  
    
    # 1 Add any persistent deformation
    for persistent_def in persistent_defs:
        d_start = datetime.strptime(str(persistent_def['def_episode_start']), '%Y%m%d')
        d_stop = datetime.strptime(str(persistent_def['def_episode_stop']), '%Y%m%d')
        r2 = Range(start = d_start, end = d_stop)
        latest_start = max(r1.start, r2.start)
        earliest_end = min(r1.end, r2.end)
        delta = (earliest_end - latest_start).days
        overlap = max(0, delta)         
        # print(f"Overlapping days: {overlap}")                                                                     # can be useful to debug/test
        if overlap > 0:
            def_predicted += ((overlap / 365.25) * persistent_def['def_rate'])                                      # convert days to years, then multiply by rate in m/year
            if persistent_def['source'] not in sources:
                sources.append(persistent_def['source'])
                        
    #print(f"{def_predicted} m")
    
    # 2: Add any transient deformation
    for transient_def in transient_defs:
        d_start = datetime.strptime(str(transient_def['def_episode_start']), '%Y%m%d')
        d_stop = datetime.strptime(str(transient_def['def_episode_stop']), '%Y%m%d')
        r2 = Range(start = d_start, end = d_stop)
        latest_start = max(r1.start, r2.start)
        earliest_end = min(r1.end, r2.end)
        delta = (earliest_end - latest_start).days
        overlap = max(0, delta)
        #print(f"Overlapping days: {overlap}")                                                                      # can be useful to debug/test
        if overlap > 0:
            def_predicted += transient_def['def_magnitude']                           # convert days to years
            if transient_def['source'] not in sources:
                sources.append(transient_def['source'])
    
    
    return def_predicted, sources

        
        
    
# def_predicted = volcnet_labeller("20180320_20190701", persistent_defs, transient_defs)
# def_predicted = volcnet_labeller("20180320_20180322", persistent_defs, transient_defs)

# def_predicted = volcnet_labeller("20141031_20210912", persistent_defs, transient_defs)
# def_predicted = volcnet_labeller("20141031_20230912", persistent_defs, transient_defs)
# def_predicted = volcnet_labeller("20101031_20230912", persistent_defs, transient_defs)

#def_predicted = volcnet_labeller("20180526_20180713", persistent_defs, transient_defs)
def_predicted, sources = volcnet_labeller("20180526_20180713", persistent_defs, transient_defs)                      # Sierra Negra co-eruptive 2018 interferograms.  
    
