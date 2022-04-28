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
import matplotlib.gridspec as gridspec

# MEG debug imports
sys.path.append("/home/matthew/university_work/python_stuff/python_scripts")
from small_plot_functions import matrix_show, pngs_to_gif, col_to_ma, quick_linegraph
import matplotlib.pyplot as plt
#plt.switch_backend('Qt5Agg')                                                               # only needed if importing this mid debug 

sys.path.append("/home/matthew/university_work/23_insar_tools")                  # 
import insar_tools
from insar_tools.plotting import xticks_every_nmonths

def remappedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
    '''
    Function to offset the median value of a colormap, and scale the
    remaining color range (i.e. truncate the colormap so that it isn't
    compressed on the shorter side) . Useful for data with a negative minimum and
    positive maximum where you want the middle of the colormap's dynamic
    range to be at zero.
    Input
    -----
      cmap : The matplotlib colormap to be altered
      start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower ofset). Should be between
          0.0 and 0.5; if your dataset mean is negative you should leave
          this at 0.0, otherwise to (vmax-abs(vmin))/(2*vmax)
      midpoint : The new center of the colormap. Defaults to
          0.5 (no shift). Should be between 0.0 and 1.0; usually the
          optimal value is abs(vmin)/(vmax+abs(vmin))
          Only got this to work with:
              1 - vmin/(vmax + abs(vmin))
      stop : Offset from highets point in the colormap's range.
          Defaults to 1.0 (no upper ofset). Should be between
          0.5 and 1.0; if your dataset mean is positive you should leave
          this at 1.0, otherwise to (abs(vmin)-vmax)/(2*abs(vmin))

      2017/??/?? | taken from stack exchange
      2017/10/11 | update so that crops shorter side of colorbar (so if data are in range [-1 100],
                   100 will be dark red, and -1 slightly blue (and not dark blue))
      '''
    import numpy as np
    import matplotlib
    import matplotlib.pyplot as plt

    if midpoint > 0.5:                                      # crop the top or bottom of the colourscale so it's not asymetric.
        stop=(0.5 + (1-midpoint))
    else:
        start=(0.5 - midpoint)


    cdict = { 'red': [], 'green': [], 'blue': [], 'alpha': []  }
    # regular index to compute the colors
    reg_index = np.hstack([np.linspace(start, 0.5, 128, endpoint=False),  np.linspace(0.5, stop, 129)])

    # shifted index to match the data
    shift_index = np.hstack([ np.linspace(0.0, midpoint, 128, endpoint=False), np.linspace(midpoint, 1.0, 129)])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)
        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))
    newcmap = matplotlib.colors.LinearSegmentedColormap(name, cdict)
    #plt.register_cmap(cmap=newcmap)
    return newcmap


#%% Settings

volcnet_dir = Path("/home/matthew/university_work/13_volcnet")

figsize = 10

#%% Get a list of Volcnet files




volcnet_files = sorted(glob.glob(str(volcnet_dir / '*.pkl')), key = os.path.getmtime)            # get the paths to the mat files from fabien

# for volcnet_file in volcnet_files[1:2]:
#     print("TESTING - only using Campi Flegrei volcnet file.  ")
# for volcnet_file in volcnet_files[11:12]:
#     print("TESTING - only using Wolf volcnet file.  ")
for volcnet_file in volcnet_files[10:11]:
    print("TESTING - only using Sierra Negra 128 volcnet file.  ")

    #print(volcnet_file)
    
    
    # 1: Open the file
    with open(volcnet_file, 'rb') as f:
        displacement_r3 = pickle.load(f)
        tbaseline_info = pickle.load(f)
        persistent_defs = pickle.load(f)
        transient_defs = pickle.load(f)

    ifg_resolution = 20
    acq_spacing = 5
   
    n_acq, ny, nx = displacement_r3['cumulative'].shape
    # if nx > ny:
    #     ifg_resolution_x = 30
    #     ifg_resolution_y = int((ny/nx) * ifg_resolution_x)
    # else:
    #     ifg_resolution_y = 30
    #     ifg_resolution_x = int((nx/ny) * ifg_resolution_y)

    
    
    cumulative_r3_small = displacement_r3['cumulative'][::acq_spacing,
                                                        np.array((np.linspace(0, ny-1, ifg_resolution)), dtype = int), :]           # downsample in time and y
    cumulative_r3_small = cumulative_r3_small[:,:, np.array((np.linspace(0, nx-1, ifg_resolution)), dtype = int)]                    # downsample in x
    

    n_acq_plot, ny_plot, nx_plot = cumulative_r3_small.shape                                                                                      # number of times in downsampled timeseries  
       
    figure = ma.zeros((n_acq_plot * ifg_resolution, n_acq_plot * ifg_resolution))                         # create a giant array of zeros to hold all the one channel images.  


    for row_n in range(n_acq_plot):                                                                         # loop down one row.  
        print(f"row {row_n}")
        cumulative_r3_small -= cumulative_r3_small[row_n,]                                                  # make the time series relative to this acquisition number.  
        for col_n in range(n_acq_plot):
            if col_n == row_n:                                                                              # if on the diagonal, make a square of zeros so it stands out
                mini_ifg = ma.array(np.zeros((ny_plot, nx_plot)), mask = np.zeros((ny_plot, nx_plot)))
            else:
                mini_ifg = cumulative_r3_small[col_n,]                                                      # for all others just get the ifg.  
            figure[row_n*ifg_resolution : (row_n + 1) *ifg_resolution,
                   col_n*ifg_resolution : (col_n + 1) *ifg_resolution,] = mini_ifg
            
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes


    fig = plt.figure(figsize=(2*figsize, figsize))
    grid = gridspec.GridSpec(12, 24, wspace=0.1, hspace=0.1)                        # divide into 2 sections, 1/5 for ifgs and 4/5 for components
    
    ax_all_ifgs = plt.Subplot(fig, grid[:12, :12])                                                # a thin but wide axes for all the thumbnail ifgs along the top to go in
    ax_1_ifg = plt.Subplot(fig, grid[:8, 12:])                                                # a thin but wide axes for all the thumbnail ifgs along the top to go in
    ax_labels = plt.Subplot(fig, grid[9:, 13:-1])                                                # a thin but wide axes for all the thumbnail ifgs along the top to go in
    
    
    cmap_mid = 1 - ma.max(figure)/(ma.max(figure) + abs(ma.min(figure)))                                     # get the ratio of the data that 0 lies at (eg if data is -15 to 5, ratio is 0.75)
    figure_cmap = remappedColorMap(plt.get_cmap('coolwarm'), start=0.0, midpoint=cmap_mid, stop=1, name='ic_colours_cent')                    # make the colours for plotting the ICs
    all_ifgs = ax_all_ifgs.imshow(figure)                                                                   # plot the giant overview of tiny interferograms.  

    def click(event):
        if event.inaxes == ax_all_ifgs:                                                                    # determine if the mouse is in the axes on the left
            if all_ifgs.contains(event):                                                               # cont is a boolean of if hoving on point, ind is a dictionary about the point being hovered over.  Note that two or more points can be in this.  
                try:
                    for ax in fig.axes[4:]:
                        ax.remove()
                except:
                    pass
                
                acq_primary = int(event.xdata / ifg_resolution) * acq_spacing                          # xdata is in pixels of the big image, figure, convert to which acquisition number this is.   
                acq_secondary = int(event.ydata / ifg_resolution) * acq_spacing                        # and in y 
                #print(f"x:{acq_primary} y:{acq_secondary}")
                ifg_array = displacement_r3['cumulative'][acq_primary] - displacement_r3['cumulative'][acq_secondary]  
                cmap_mid = 1 - ma.max(ifg_array)/(ma.max(ifg_array) + abs(ma.min(ifg_array)))                                     # get the ratio of the data that 0 lies at (eg if data is -15 to 5, ratio is 0.75)
                colours_cent = remappedColorMap(plt.get_cmap('coolwarm'), start=0.0, midpoint=cmap_mid, stop=1, name='ic_colours_cent')                    # make the colours for plotting the ICs
                ifg = ax_1_ifg.imshow(ifg_array, cmap = colours_cent)   # draw that ifg.  
                
                primary = datetime.strptime(tbaseline_info['acq_dates'][acq_primary], '%Y%m%d')
                secondary = datetime.strptime(tbaseline_info['acq_dates'][acq_secondary], '%Y%m%d')
                tbaseline = (primary - secondary).days
                
                ax_1_ifg.set_title(f"{tbaseline_info['acq_dates'][acq_secondary]}_{tbaseline_info['acq_dates'][acq_primary]} ({tbaseline} days)")
                
                
                # axins2 = inset_axes(ax_1_ifg, width="100%", height="100%",                                                             # what fraction of the bounding box to take up
                #                     bbox_to_anchor=(0.25, -0.2, 0.5, 0.07), bbox_transform=ax_1_ifg.transAxes)                         # x start y start x width y height    
                axins2 = inset_axes(ax_1_ifg, width="100%", height="100%", loc = 'upper left',                                                             # what fraction of the bounding box to take up
                                    bbox_to_anchor=(0.05, 0.95, 0.3, 0.05), bbox_transform=ax_1_ifg.transAxes)                         # x start y start x width y height    
                cb = fig.colorbar(ifg, cax = axins2, orientation = 'horizontal')
                cb.set_label("LOS displacement (m)")
                
                # update the tick labels to be lon lat
                xtick_labels = []
                for tick in ax_1_ifg.get_xticks():
                    try:                                                                                                    # ticks can extend past data 
                        xtick_labels.append(np.round(displacement_r3['lons'][-1, int(tick)], 2))        
                    except:
                        xtick_labels.append('')
                ax_1_ifg.xaxis.set_ticklabels(xtick_labels, rotation = 315, ha = 'left')
                
                ytick_labels = []
                for tick in ax_1_ifg.get_yticks():
                    try:                                                                                                    # ticks can extend past data 
                        ytick_labels.append(np.round(displacement_r3['lats'][int(tick), 0], 2))        
                    except:
                        ytick_labels.append('')
                ax_1_ifg.yaxis.set_ticklabels(ytick_labels)
                
                
                
            else:                                                                       # else not on a point
                pass
        else:                                                                           # else not in the axes
            pass
        fig.canvas.draw_idle()
    
    

    fig.canvas.mpl_connect("button_press_event", click)                                # connect the figure and the function, note that done when the mouse clicks.  
    ax_all_ifgs.set_xlabel('Primary acquisition')
    ax_all_ifgs.set_ylabel('Secondary acquisition')
    # Set the tick labels to be in dates and not pixels of the giant figure.  
    xtick_labels = []
    for tick in ax_all_ifgs.get_xticks():
        try:                                                                                                    # ticks can extend past data 
            xtick_labels.append(tbaseline_info['acq_dates'][int((tick / ifg_resolution) * acq_spacing)])        
        except:
            xtick_labels.append('')
    ax_all_ifgs.xaxis.set_ticklabels(xtick_labels, rotation = 315, ha = 'left')
    ax_all_ifgs.xaxis.tick_top()
    ax_all_ifgs.xaxis.set_label_position('top')
    ytick_labels = []
    for tick in ax_all_ifgs.get_yticks():
        try:                                                                                                    # ticks can extend past data 
            ytick_labels.append(tbaseline_info['acq_dates'][int((tick / ifg_resolution) * acq_spacing)])        
        except:
            ytick_labels.append('')
    ax_all_ifgs.yaxis.set_ticklabels(ytick_labels)
    
    fig.add_subplot(ax_all_ifgs)                                                                   # add to figure
    fig.add_subplot(ax_1_ifg)                                                                   # add to figure
    
    
    #Add the label information at the bottom.  
    
    fig.add_subplot(ax_labels)                                                                   # add to figure
    
    # ax_labels = inset_axes(ax_all_ifgs, width="100%", height="100%",                                                             # what fraction of the bounding box to take up
    #                        bbox_to_anchor=(0., -0.1, 1., 0.1), bbox_transform=ax_all_ifgs.transAxes)                         # x start y start x width y height    
    ax_transient = ax_labels.twinx()
    
    d0_dt = datetime.strptime(tbaseline_info['acq_dates'][0], '%Y%m%d')
    ax_labels.set_xlim(tbaseline_info['baselines_cumulative'][0], tbaseline_info['baselines_cumulative'][-1])
    
    for persistent_def in persistent_defs:
        colour = 'tab:orange'
        x_start = (datetime.strptime(str(persistent_def["def_episode_start"]), '%Y%m%d') - d0_dt).days
        x_stop = (datetime.strptime(str(persistent_def["def_episode_stop"]), '%Y%m%d') - d0_dt).days
        ax_transient.plot((x_start, x_stop), (persistent_def["def_rate"], persistent_def["def_rate"]), c = colour)
        
    for transient_def in transient_defs:
        colour = 'tab:blue'
        x_start = (datetime.strptime(str(transient_def["def_episode_start"]), '%Y%m%d') - d0_dt).days
        x_stop = (datetime.strptime(str(transient_def["def_episode_stop"]), '%Y%m%d') - d0_dt).days
        ax_transient.plot((x_start, x_stop), (transient_def["def_magnitude"], transient_def["def_magnitude"]), c = colour)



    xticks_every_nmonths(ax_labels, tbaseline_info['acq_dates'][0], tbaseline_info['baselines_cumulative'], include_tick_labels = True,                  # update x ticks, but with labels.  
                         major_ticks_n_months = 12, minor_ticks_n_months = 3)

    ax_transient.set_ylabel('Transient\ndeformation (m)', color=colour, fontsize= 10)  # we already handled the x-label with ax1
    ax_transient.tick_params(axis='y', labelcolor=colour)
    ax_labels.set_ylabel('Persistent\ndeformation (m/yr)', fontsize = 10)
        # pdb.set_trace()
        
    
    #fig.add_subplot(ax_labels)                                                                   # add to figure
    
    
#%% create labels

    
# make an ifg.  
# check if transient

# check if persistent def exists.  
# estimate rate from persisent.  
# compare to threshold
# assign label.  

def volcnet_labeller(ifg_name, persistent_defs, transient_defs):
    """
    """
    from collections import namedtuple
    
    Range = namedtuple('Range', ['start', 'end'])
    
    acq_start_dt = datetime.strptime(ifg_name[:8], '%Y%m%d')
    acq_stop_dt = datetime.strptime(ifg_name[9:], '%Y%m%d')
    tbaseline = (acq_stop_dt - acq_start_dt).days
    r1 = Range(start = acq_start_dt, end = acq_stop_dt)
    def_predicted = 0.
    
    # 1 Add any persistent deformation
    for persistent_def in persistent_defs:
        d_start = datetime.strptime(str(persistent_def['def_episode_start']), '%Y%m%d')
        d_stop = datetime.strptime(str(persistent_def['def_episode_stop']), '%Y%m%d')
        r2 = Range(start = d_start, end = d_stop)
        latest_start = max(r1.start, r2.start)
        earliest_end = min(r1.end, r2.end)
        delta = (earliest_end - latest_start).days
        overlap = max(0, delta)
        # print(f"Overlapping days: {overlap}")
        def_predicted += ((overlap / 365.25) * persistent_def['def_rate'])                           # convert days to years
                        
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
        print(f"Overlapping days: {overlap}")
        if overlap > 0:
            def_predicted += transient_def['def_magnitude']                           # convert days to years
    
    print(def_predicted)
    return def_predicted

        
        
    
# def_predicted = volcnet_labeller("20180320_20190701", persistent_defs, transient_defs)
# def_predicted = volcnet_labeller("20180320_20180322", persistent_defs, transient_defs)

# def_predicted = volcnet_labeller("20141031_20210912", persistent_defs, transient_defs)
# def_predicted = volcnet_labeller("20141031_20230912", persistent_defs, transient_defs)
# def_predicted = volcnet_labeller("20101031_20230912", persistent_defs, transient_defs)

#def_predicted = volcnet_labeller("20180526_20180713", persistent_defs, transient_defs)
def_predicted = volcnet_labeller("20180526_20180713", persistent_defs, transient_defs)


