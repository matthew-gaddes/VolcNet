#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 10:58:58 2022

@author: matthew


################################################################################################################################

This function builds the VolcNet database from various input data (e.g. Sentinel-1 time series, COMET LiCS time series etc.)

Contains:
    
    Lin Shen                            XXXXXXXX
    Matthew Gaddes                      LiCSAR + LiCSBAS time series for various volcanoes.  
    Marco Bagnardi + Matthew Gaddes     Galapagos 
    Fabien Albino                       Agung.  

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


#ICASAR_path = Path("/home/matthew/university_work/15_my_software_releases/ICASAR-2.11.0/")                               # 
ICASAR_path = Path('/home/matthew/university_work/01_blind_signal_separation_python/13_ICASAR/ICASAR_GitHub')          # development
#print(f"USING THE GITHUB VERSION OF ICASAR")
sys.path.append(str(ICASAR_path))
import icasar
from icasar.icasar_funcs import LiCSBAS_to_ICASAR
from icasar.aux import r2_arrays_to_googleEarth


sys.path.append("/home/matthew/university_work/23_insar_tools")                  # 
import insar_tools
from insar_tools.general import r2_to_r3, reference_r3_ts
from insar_tools.temporal_baselines import acquisitions_from_ifg_dates, daisy_chain_from_acquisitions, baselines_from_names
from insar_tools.open_data import open_fabien_cf_data


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
        list_of_strings = config.get(deformation_event, 'def_polygon').strip('][').split('),')                                      # bit of a horrible string that looks like a list of tuples, convert to list of messy strings
        event['def_polygon'] = []                                                                                                   # initiate
        for list_of_string in list_of_strings:                                                                                      # loop through
            event['def_polygon'].append((float(list_of_string.strip(' ()').split(', ')[0]), float(list_of_string.strip(' ()').split(', ')[1])))       # convert string that looks like tuple to a tuple.  
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

#############################




#%% Things to set

licsbas_ts_dir = Path("/home/matthew/university_work/data/00_LiCSBAS_time_series")                              # directory that LiCSBAS time series are stored in.  
licsbas_labels_dir = Path("/home/matthew/university_work/13_volcnet/raw_annotation_data/licsbas")
mask_type = 'licsbas'                                                                                           # can be 'dem' (more pixels) or 'licsbas'

marco_galapagos_dir = Path("/home/matthew/university_work/data/01_galapagos/pre_licsar_timeseries")
marco_galapagos_labels_dir = Path("/home/matthew/university_work/13_volcnet/raw_annotation_data/marco_galapagos")
s1_lambda = 0.0546573                                                                                               # wavelength, in metres, needed to convert galapagos data from rads to m

fabien_agung_dir = Path("/home/matthew/university_work/data/fabien_agung")
fabien_agung_labels_dir = Path("/home/matthew/university_work/13_volcnet/raw_annotation_data/fabien_agung")

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


#%% Time series that were processed by LiCSBAS.  

for licsbas_ts in licsbas_tss:
    
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
    #r2_arrays_to_googleEarth(displacement_r3['cumulative'][-1,][np.newaxis,], displacement_r3['lons'], displacement_r3['lats'], 'IC', out_folder = './', kmz_filename = 'last_cumulative_ifg')                              # note that lons and lats should be rank 2 (ie an entry for each pixel in the ifgs)
    
    # open the aux data
    persistent_defs, transient_defs = read_volcnet_label(licsbas_labels_dir / (licsbas_ts + '.txt'))
    
    # remove some attributes as we don't have them for the Marco Galapagos data and want to be consistent.  
    del displacement_r3['E'], displacement_r3['N'], displacement_r3['U']
    
    # possibly rename from MEG licsbas time series names to Volcnet names (frame + volcname)
    if licsbas_ts == "002A_05136_020502_azores_crop":
        licsbas_ts = "002A_05136_020502_azores_crop_sao_jorge"
    elif licsbas_ts == '022D_04826_121209':
        licsbas_ts = "022D_04826_121209_campi_flegrei"
    elif licsbas_ts == '022D_04826_121209_vesuvius_crop':
        licsbas_ts = "022D_04826_121209_vesuvius"
    elif licsbas_ts == '082D_05128_030500_azores_crop':
        licsbas_ts = "082D_05128_030500_azores_sao_jorge"
    elif licsbas_ts == '083D_12636_131313_domuyo_rationalized_v2':
        licsbas_ts = '083D_12636_131313_domuyo'
    elif licsbas_ts == '124D_04854_171313_licsbas_example_extended':
        licsbas_ts = '124D_04854_171313_campi_flegrei'
    elif licsbas_ts == '169D_00001_020800_rationalized':
        licsbas_ts = "169D_00001_020800_la_plama"
        
    
    
    # write to a file
    with open(f"{licsbas_ts}.pkl", 'wb') as f:                                                                     # Or open the products from a previous ICASAR run.  
        pickle.dump(displacement_r3, f)
        pickle.dump(tbaseline_info, f)
        pickle.dump(persistent_defs, f)
        pickle.dump(transient_defs, f)
    f.close()                                                                                                                           
    
    
#%% Marco Bagnardi Galapgos time series.  

ref_areas = {"106A_09090_000404_sierra_negra"   : (300, 50, 50, 50),                  # x start, y start, x window size, y window size
              "106A_09090_000404_cerro_azul"     : (240, 35, 50, 50),
              "128D_09016_110500_sierra_negra"   : (260, 31, 50, 50),
              "128D_09016_110500_wolf"           : (337, 390, 50, 50),
              "128D_09016_110500_cerro_azul"     : (270, 32, 50, 50),
              "106A_09090_000404_wolf"           : (337, 390, 50, 50)}


# get a list of the files (but not the lon lat ones)
galapagos_files = glob.glob(str(marco_galapagos_dir / "*"))                                                # get the paths to each template file
ll_file_indexes = []
for file_n, galapagos_file in enumerate(galapagos_files):
    if galapagos_file[-7:-4] == '_ll':
        ll_file_indexes.append(file_n)
for ll_file_index in sorted(ll_file_indexes)[::-1]:
    del galapagos_files[ll_file_index]
if len(galapagos_files) == 0:
    raise Exception("No galapagos files were found.  Unable to proceed.")
  
#pdb.set_trace()
  
for galapagos_file in galapagos_files[4:5]:
    
    ts_name = galapagos_file.split('/')[-1].split('.')[0]
    
    print(f"Opening the LiCSAR ifgs processed by Marco Bagnardi and Matthew Gaddes for {ts_name}...")
    # 0: Open the data:
    with open(galapagos_file, 'rb') as f:
        mask_water = pickle.load(f)    
        mask_coh_water = pickle.load(f)    
        ifg_names = pickle.load(f)    
        phUnw_unreferenced = pickle.load(f)                                                          # still in rads and positive is range increase so subsidence.  
        dem = pickle.load(f)        
    
    with open(galapagos_file[:-4] + "_ll.pkl", 'rb') as f:
        lons = pickle.load(f)
        lats = pickle.load(f)
    
    
    # 0 Reference time series
    
    phUnw = reference_r3_ts(phUnw_unreferenced, ref_areas[ts_name][0], ref_areas[ts_name][1], ref_areas[ts_name][2], ref_areas[ts_name][3])
        
    # 1 convert format to Volcnet Style
    displacement_r3 = {'lons' : lons,
                        'lats' : lats,
                        'dem'  : ma.masked_array(dem, mask = mask_coh_water),
                        'mask' : mask_coh_water}
    
    # 1a: get only the daisy chain ifgs.  
    acq_dates = acquisitions_from_ifg_dates(ifg_names)                                                                               # get the S1 acquisitions
    daisy_chains = daisy_chain_from_acquisitions(acq_dates)                                                                          # make a simple daisy chain time series.  
    cumulative_r3 = np.zeros(((1, phUnw.shape[1], phUnw.shape[2])))                                                                  # initiate to store the daisy chain
    for daisy_chain in daisy_chains:                                                                                                 # loop through each one
        ifg_n = ifg_names.index(daisy_chain)                                                                                         # find which ifg_n it is.  
        cumulative_r3 = np.concatenate((cumulative_r3, phUnw[ifg_n:ifg_n+1,]), axis = 0)                                             # and add it to the big array.  
    cumulative_r3 = np.cumsum(cumulative_r3, axis = 0)                                                                               # make cumulative
    cumulative_r3 *= (s1_lambda / (4*np.pi))                                                                                         # convert from rads to m                
    cumulative_r3 *= (-1)                                                                                                            # positive is range increase so subsidence.  for m, convert to los displacment (up is positive)
    cumulative_r3_ma = ma.masked_array(cumulative_r3, np.repeat(mask_coh_water[np.newaxis,], cumulative_r3.shape[0], axis = 0))      # mask water and incoherence 
    displacement_r3['cumulative'] = cumulative_r3_ma
    incremental_r3 = ma.diff(cumulative_r3_ma, axis = 0)
    
    # 2: get the time info
    tbaseline_info = {}
    tbaseline_info['acq_dates'] = acq_dates
    baselines_dc = baselines_from_names(daisy_chains)
    baselines_cs = np.cumsum(baselines_dc)
    baselines_cs = np.concatenate((np.array([0]), baselines_cs))
    tbaseline_info['baselines_cumulative'] = baselines_cs
    
    ############################## optional to inspect data
    if False:
        # Possibly plot a time series
#        quick_linegraph([cumulative_r3[:, 200, 200], cumulative_r3[:, 310, 50]], xvals = baselines_cs)                                                         # 106A sierra negra
        #quick_linegraph([cumulative_r3[:, 304, 276], cumulative_r3[:, ref_areas[ts_name][1], ref_areas[ts_name][0]]], xvals = baselines_cs)                    # 106A cero azul
        #quick_linegraph([cumulative_r3[:, 211, 218], cumulative_r3[:, ref_areas[ts_name][1], ref_areas[ts_name][0]]], xvals = baselines_cs)                     # 128D Sierra Negra
        quick_linegraph([cumulative_r3[:, 295, 341], cumulative_r3[:, ref_areas[ts_name][1], ref_areas[ts_name][0]]])#, xvals = baselines_cs)                     # 128D cerro azul
        t_start = 935           
        t_stop = 1300
        y_start = 1.03
        y_stop = 2.5
        gradient = (y_stop - y_start) / ((t_stop - t_start)/365)
        
        
        # # to inspect ifgs.  
        for i in np.arange(10, 20):
        #for i in np.arange(1, len(acq_dates)):
            f, axes = plt.subplots(1,2, figsize = (12,6))
            matrix_show(cumulative_r3_ma[i,], title = f"cum ifg {i}: {acq_dates[0]} - {acq_dates[i]}", fig = f, ax= axes[0])
            matrix_show(incremental_r3[i-1,], title = f"inc ifg {i-1}: {daisy_chains[i-1]}", fig = f, ax= axes[1])
        
        
        # # possibly output to google earth to check
        r2_arrays_to_googleEarth(displacement_r3['cumulative'][-1,][np.newaxis,], displacement_r3['lons'], displacement_r3['lats'], 'IC', out_folder = './', kmz_filename = 'last_cumulative_ifg')                              # note that lons and lats should be rank 2 (ie an entry for each pixel in the ifgs)
        r2_arrays_to_googleEarth(incremental_r3[8:11,], displacement_r3['lons'], displacement_r3['lats'], 'IC', out_folder = './', kmz_filename = 'last_cumulative_ifg')                              # note that lons and lats should be rank 2 (ie an entry for each pixel in the ifgs)
    ############################## end optional to inspect data
  

    # open the aux data
    persistent_defs, transient_defs = read_volcnet_label(marco_galapagos_labels_dir / (ts_name + '.txt'))

    
    
    # 2: write to a file
    with open(f"{ts_name}.pkl", 'wb') as f:                                                                     # Or open the products from a previous ICASAR run.  
        pickle.dump(displacement_r3, f)
        pickle.dump(tbaseline_info, f)
        pickle.dump(persistent_defs, f)
        pickle.dump(transient_defs, f)
    f.close()                                                                                                                           

    



#%% Fabien Albino Agung data.  


def load_fabien_defo_file(defo_file, ifg_dates_file):
    """A function to load one of Fabien's tiem series files.  Dates of the incremental interferograms 
    are in a text file, and the other info is in .mat file.  
    He treats icoherent pixels the same as water, and sets these to 0?  This function uses that create a mask
    that removes both types of pixel.  
    Inputs:
        defo_file | .mat file | file from Fabien.  
        ifg_dates_file | .txt file | ifg dates, one on each line.  Of form YYYYMMDD_YYYYMMDD
    """
    import scipy.io as sio
    
    # 1: open the text file of dates
    print(f'Defo file: {defo_file}, Dates file: {ifg_dates_file}')
    f = open(ifg_dates_file, "r")
    all_dates = f.read()
    all_dates = all_dates.split('\n')
    ifg_names = []
    for ifg_n in np.arange(1, len(all_dates)):
        ifg_names.append(f'{all_dates[ifg_n-1]}_{all_dates[ifg_n]}')
    
    
    # 2: open the data
    mat_contents = sio.loadmat(defo_file)
    incremental = mat_contents['incremental']                                                                   # these are a rank 3 of data (heigh x wdith x n_ifgs)
    lons = mat_contents['LON']                                                                                  # rank 1 (ie row vector?)
    lats = mat_contents['LAT']                                                                                  # rank 1 (ie a column vector?)    
    

    # start processing the data
    incremental = np.rollaxis(incremental, 2, 0)                                                                        # get time as first dimension
    mask = np.where(incremental[0,:,:] == 0, np.ones(incremental[0,:,:].shape), np.zeros(incremental[0,:,:].shape))
    print(f'One pixel is {111000 * (lats[0,0] - lats[0,1])} metres')
    
    # make a masked array
    n_ifgs = incremental.shape[0]
    phUnw = ma.zeros(incremental.shape)                      # rank 3 array, in the Keras style
    for ifg_n in range(n_ifgs):
        phUnw[ifg_n,:,:] = ma.array(incremental[ifg_n,:,:], mask = mask)
    
    # check some size
    print(f'There are {incremental.shape[0]} ifgs, and {len(ifg_names)} ifg names')
    return phUnw, mask, ifg_names, lons, lats





# 1: Agung descedingls -lh
# deformation isn't very clear

# 2: Agung ascneding

phUnw_ma, mask_water_coh, ifg_names, lons_r1, lats_r1 = load_fabien_defo_file(fabien_agung_dir / "ag_txxxA_combined.mat",
                                                                              fabien_agung_dir / "ascending_date.txt")  

lons, lats = np.meshgrid(lons_r1, lats_r1)

displacement_r3 = {'lons' : lons,
                   'lats' : lats,
                   'mask' : mask_water_coh}

# time info
acq_dates = acquisitions_from_ifg_dates(ifg_names)                                                                               # get the S1 acquisitions
daisy_chains = daisy_chain_from_acquisitions(acq_dates)                                                                          # make a simple daisy chain time series.  
baselines_dc = baselines_from_names(daisy_chains)
baselines_cs = np.cumsum(baselines_dc)
baselines_cs = np.concatenate((np.array([0]), baselines_cs))

tbaseline_info = {'acq_dates' : acq_dates,
                  'baselines_cumulative' : baselines_cs}

# LOS displacment info
cumulative_r3 = ma.concatenate((np.zeros((1, phUnw_ma.shape[1], phUnw_ma.shape[2])), np.cumsum(phUnw_ma, axis = 0)), axis = 0)
incremental_r3 = ma.diff(cumulative_r3, axis = 0)

displacement_r3['cumulative'] = cumulative_r3
############################## optional to inspect data
if False:
    # Possibly plot a time series
  
    quick_linegraph([cumulative_r3[:, 337, 643], cumulative_r3[:, 295, 341] ])#, xvals = baselines_cs)                     # 128D cerro azul
    quick_linegraph([cumulative_r3[:, 337, 647], cumulative_r3[:, 295, 341] ], xvals = baselines_cs)                     # 128D cerro azul
    t_start = 935           
    t_stop = 1300
    y_start = 1.03
    y_stop = 2.5
    gradient = (y_stop - y_start) / ((t_stop - t_start)/365)
    
    
    # # to inspect ifgs.  
    for i in np.arange(8, 15):
    # for i in np.arange(1, len(acq_dates)):
        f, axes = plt.subplots(1,2, figsize = (12,6))
        matrix_show(cumulative_r3[i,], title = f"cum ifg {i}: {acq_dates[0]} - {acq_dates[i]}", fig = f, ax= axes[0])
        matrix_show(incremental_r3[i-1,], title = f"inc ifg {i-1}: {daisy_chains[i-1]}", fig = f, ax= axes[1])
    
    
    # # possibly output to google earth to check
    r2_arrays_to_googleEarth(displacement_r3['cumulative'][-1,][np.newaxis,], displacement_r3['lons'], displacement_r3['lats'], 'IC', out_folder = './', kmz_filename = 'last_cumulative_ifg')                              # note that lons and lats should be rank 2 (ie an entry for each pixel in the ifgs)
    r2_arrays_to_googleEarth(incremental_r3[8:11,], displacement_r3['lons'], displacement_r3['lats'], 'IC', out_folder = './', kmz_filename = 'last_cumulative_ifg')                              # note that lons and lats should be rank 2 (ie an entry for each pixel in the ifgs)
############################## end optional to inspect data
  

# open the aux data
persistent_defs, transient_defs = read_volcnet_label(fabien_agung_labels_dir / "156A_09814_081406_agung.txt")    


# 2: write to a file
with open(f"156A_09814_081406_agung.pkl", 'wb') as f:                                                                     # Or open the products from a previous ICASAR run.  
    pickle.dump(displacement_r3, f)
    pickle.dump(tbaseline_info, f)
    pickle.dump(persistent_defs, f)
    pickle.dump(transient_defs, f)
f.close()                                                                                                                               