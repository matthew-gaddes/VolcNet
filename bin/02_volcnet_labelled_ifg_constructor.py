#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 15:10:50 2022

@author: matthew

################################################################################################################################

This function creates labelled interferograms from the VolcNet data.  


################################################################################################################################

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
from volcnet.plotting import volcnet_ts_visualiser, plot_volcnet_files_labels
from volcnet.labelling import label_volcnet_ifg, label_volcnet_files


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
volcnet_def_min = 0.05                                                              # in m, less than this is just considered as containing noise.  

figsize = 10

#%% Visualise each time series (can be slow!)



volcnet_files = sorted(glob.glob(str(volcnet_dir / '*.pkl')))            # get the paths to the mat files from fabien

# for volcnet_file in volcnet_files[1:2]:
#     print("TESTING - only using Campi Flegrei volcnet file.  ")
# for volcnet_file in volcnet_files[11:12]:
#     print("TESTING - only using Wolf volcnet file.  ")
for volcnet_file in volcnet_files[11:12]:
    print("TESTING - only using Sierra Negra 128 volcnet file.  ")
# for volcnet_file in volcnet_files:
    
    
    # 1: Open the file
    with open(volcnet_file, 'rb') as f:
        displacement_r3 = pickle.load(f)
        tbaseline_info = pickle.load(f)
        persistent_defs = pickle.load(f)
        transient_defs = pickle.load(f)

    print(volcnet_file)
    volcnet_ts_visualiser(displacement_r3, tbaseline_info, persistent_defs, transient_defs, acq_spacing = 1, ifg_resolution = 20, figsize_height = 10,
                          labelling_function = label_volcnet_ifg, title = Path(volcnet_file).parts[-1].split('.')[0])                                               # last part gets the frane name from the path, regardless of operating system.  
    
    

#%% Visualise the whole database

labels_dyke, labels_sill, labels_atmo = label_volcnet_files(volcnet_files, def_min = volcnet_def_min)               # label all the VolcNet data - slow! 
# plot_volcnet_files_labels(volcnet_files, labels_dyke, labels_sill, labels_atmo)                                     # figure with two subplots showing how many of each label etc.  




# with open('fig_6_data.pkl', 'wb') as f:                                       # Used in CNN paper to show breakdown of labels by type and volcano etc.  
#     pickle.dump(volcnet_files, f)
#     pickle.dump(labels_dyke, f)
#     pickle.dump(labels_sill, f)
#     pickle.dump(labels_atmo, f)
    
#%% Test labelling of a single interferogram.  

        
    
# def_predicted = volcnet_labeller("20180320_20190701", persistent_defs, transient_defs)
# def_predicted = volcnet_labeller("20180320_20180322", persistent_defs, transient_defs)

# def_predicted = volcnet_labeller("20141031_20210912", persistent_defs, transient_defs)
# def_predicted = volcnet_labeller("20141031_20230912", persistent_defs, transient_defs)
# def_predicted = volcnet_labeller("20101031_20230912", persistent_defs, transient_defs)

#def_predicted = volcnet_labeller("20180526_20180713", persistent_defs, transient_defs)
#def_predicted, sources, def_location = volcnet_labeller("20180526_20180713", persistent_defs, transient_defs)                      # Sierra Negra co-eruptive 2018 interferograms.  
#def_predicted, sources, def_location = label_volcnet_ifg("20180713_20180526", persistent_defs, transient_defs)                      # Sierra Negra co-eruptive 2018 interferograms.  


