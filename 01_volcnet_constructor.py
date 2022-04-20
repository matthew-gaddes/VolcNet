#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 10:58:58 2022

@author: matthew

Thoughts:
    
- Loop through each LiCSBAS time series. 
    - Open LiCSBAS data.  
    - Open aux data.  
        - location
        - persistent rate for deformation
        - isolated time for deformation.  
    - save as a hdf file?  

Also Fabien Agung
And Marco Sierra Negra?

"""

import numpy as np
import matplotlib.pyplot as plt
import sys
from pathlib import Path
from copy import deepcopy
import pdb
import pickle


#ICASAR_path = Path("/home/matthew/university_work/15_my_software_releases/ICASAR-2.11.0/")                               # 
ICASAR_path = Path('/home/matthew/university_work/01_blind_signal_separation_python/13_ICASAR/ICASAR_GitHub')          # development
#print(f"USING THE GITHUB VERSION OF ICASAR")
sys.path.append(str(ICASAR_path))
import icasar
from icasar.icasar_funcs import LiCSBAS_to_ICASAR
from icasar.aux import r2_arrays_to_googleEarth


sys.path.append("/home/matthew/university_work/23_insar_tools")                  # 
import insar_tools
from insar_tools.general import r2_to_r3


# MEG debug imports
sys.path.append("/home/matthew/university_work/python_stuff/python_scripts")
from small_plot_functions import matrix_show, pngs_to_gif, col_to_ma, quick_linegraph
import matplotlib.pyplot as plt
#plt.switch_backend('Qt5Agg')                                                               # only needed if importing this mid debug 



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
        
        if deformation_event[:10] == 'persistent':                                                          # then open only persistent event paramters
            event['def_rate'] = float(config.get(deformation_event, 'def_rate'))       # 
            persistent_deformation_events.append(event)
            
        else:                                                                                               # or transient event parameters.  
            event['def_magnitude'] = float(config.get(deformation_event, 'def_magnitude'))       # 
            transient_deformation_events.append(event)
        
    return persistent_deformation_events, transient_deformation_events


#%% Things to set

licsbas_ts_dir = Path("/home/matthew/university_work/data/00_LiCSBAS_time_series")                              # directory that LiCSBAS time series are stored in.  
licsbas_labels_dir = Path("/home/matthew/university_work/13_volcnet/raw_annotation_data/licsbas")
mask_type = 'licsbas'                                                                                           # can be 'dem' (more pixels) or 'licsbas'

#%% Define Licsbas volcanoes that are included:
    
licsbas_tss = [#'002A_05136_020502_azores',
              '002A_05136_020502_azores_crop',                          # no deformation visible, poor coherence
              #'018A_12668_131313_domuyo',                           # doesn't have any data in
              '022D_04826_121209',                                      # campi flegrei
              #'022D_04826_121209_vesuvius',
              '022D_04826_121209_vesuvius_crop',
              #'060A_00001_030604_rationalized',
              '079D_07694_131313_erta_ale',                         # erption in, info taken from Chris Moore paper.  
              #'082D_05128_030500_azores',
              '082D_05128_030500_azores_crop',                      # no deformation visible, poor coherence.  
              #'083D_12636_131313_domuyo',                                   # appears to be empty
              '083D_12636_131313_domuyo_rationalized_v2',
              '124D_04854_171313_licsbas_example_extended',                     # campi flegrei
              #'128D_09016_110500_sierra_negra',                                     # no data in cum.h5
              #'128D_09016_110500_sierra_negra_bright',
              '169D_00001_020800_rationalized']                                 # la palma


#%% 

for licsbas_ts in licsbas_tss:
    #print("TESTING - CAMPI FLEGREI ONLY")
    
    # open the LiCSBAS data.  
    print(f"Opening the LiCSBAS results for {licsbas_ts}...")
    displacement_r2, tbaseline_info, ref_xy = LiCSBAS_to_ICASAR(licsbas_ts_dir / licsbas_ts,
                                                                figures=True, ref_area = True, mask_type = mask_type)    # open the h5 file produced by LiCSBAS, lons and lats are in geocode info and same resolution as the ifgs
    
    del tbaseline_info['ifg_dates'], tbaseline_info['baselines'], displacement_r2['incremental']                        # remove some data formats that are not needed (as simples to be consistent and keep only the cumulative information)
    
    # convert it to rank 3 (ie time lat lon)
    displacement_r3 = deepcopy(displacement_r2)                                                                         # make a dict for the rank 3 data
    del displacement_r3['cumulative']                                                                                   # remove this as it's the only one that's still rank 2 (ie row vectors)    
    displacement_r3['cumulative'] = r2_to_r3(displacement_r2['cumulative'], displacement_r2['mask'])                     # and remake as rank 3 (ie times x ny x nx)
    
    # possibly output to google earth to check
    #r2_arrays_to_googleEarth(displacement_r3['cumulative'][-1,][np.newaxis,], displacement_r2['lons'], displacement_r2['lats'], 'IC', out_folder = './', kmz_filename = 'last_cumulative_ifg')                              # note that lons and lats should be rank 2 (ie an entry for each pixel in the ifgs)
    
    # open the aux data
    persistent_def, transient_def = read_volcnet_label(licsbas_labels_dir / (licsbas_ts + '.txt'))
    
    # write to a file
    with open(f"{licsbas_ts}.pkl", 'wb') as f:                                                                     # Or open the products from a previous ICASAR run.  
        pickle.dump(displacement_r3, f)
        pickle.dump(tbaseline_info, f)
        pickle.dump(persistent_def, f)
        pickle.dump(transient_def, f)
    f.close()                                                                                                                           
    
    
    
    
    
    
    