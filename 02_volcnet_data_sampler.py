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

# MEG debug imports
sys.path.append("/home/matthew/university_work/python_stuff/python_scripts")
from small_plot_functions import matrix_show, pngs_to_gif, col_to_ma, quick_linegraph
import matplotlib.pyplot as plt
#plt.switch_backend('Qt5Agg')                                                               # only needed if importing this mid debug 


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

for volcnet_file in volcnet_files[1:2]:
    print("TESTING - only using Campi Flegrei volcnet file.  ")
    print(volcnet_file)
    
    
    # 1: Open the file
    with open(volcnet_file, 'rb') as f:
        displacement_r3 = pickle.load(f)
        tbaseline_info = pickle.load(f)
        persistent_def = pickle.load(f)
        transient_def = pickle.load(f)

    ifg_resolution = 10
    acq_spacing = 5
   
    n_acq, ny, nx = displacement_r3['cumulative'].shape
    # if nx > ny:
    #     ifg_resolution_x = 30
    #     ifg_resolution_y = int((ny/nx) * ifg_resolution_x)
    # else:
    #     ifg_resolution_y = 30
    #     ifg_resolution_x = int((nx/ny) * ifg_resolution_y)

    from datetime import datetime
    
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
            



    fig, axes = plt.subplots(1,2, figsize=(2*figsize, figsize))
    cmap_mid = 1 - ma.max(figure)/(ma.max(figure) + abs(ma.min(figure)))                                     # get the ratio of the data that 0 lies at (eg if data is -15 to 5, ratio is 0.75)
    figure_cmap = remappedColorMap(plt.get_cmap('coolwarm'), start=0.0, midpoint=cmap_mid, stop=1, name='ic_colours_cent')                    # make the colours for plotting the ICs
    all_ifgs = axes[0].imshow(figure)                                                                   # plot the giant overview of tiny interferograms.  

    def click(event):
        if event.inaxes == axes[0]:                                                                    # determine if the mouse is in the axes on the left
            if all_ifgs.contains(event):                                                               # cont is a boolean of if hoving on point, ind is a dictionary about the point being hovered over.  Note that two or more points can be in this.  
                try:
                    for ax in fig.axes[2:]:
                        ax.remove()
                except:
                    pass
                
                acq_primary = int(event.xdata / ifg_resolution) * acq_spacing                          # xdata is in pixels of the big image, figure, convert to which acquisition number this is.   
                acq_secondary = int(event.ydata / ifg_resolution) * acq_spacing                        # and in y 
                #print(f"x:{acq_primary} y:{acq_secondary}")
                ifg_array = displacement_r3['cumulative'][acq_primary] - displacement_r3['cumulative'][acq_secondary]  
                cmap_mid = 1 - ma.max(ifg_array)/(ma.max(ifg_array) + abs(ma.min(ifg_array)))                                     # get the ratio of the data that 0 lies at (eg if data is -15 to 5, ratio is 0.75)
                colours_cent = remappedColorMap(plt.get_cmap('coolwarm'), start=0.0, midpoint=cmap_mid, stop=1, name='ic_colours_cent')                    # make the colours for plotting the ICs
                ifg = axes[1].imshow(ifg_array, cmap = colours_cent)   # draw that ifg.  
                
                primary = datetime.strptime(tbaseline_info['acq_dates'][acq_primary], '%Y%m%d')
                secondary = datetime.strptime(tbaseline_info['acq_dates'][acq_secondary], '%Y%m%d')
                tbaseline = (primary - secondary).days
                
                axes[1].set_title(f"{tbaseline_info['acq_dates'][acq_secondary]}_{tbaseline_info['acq_dates'][acq_primary]} ({tbaseline} days)")
                
                from mpl_toolkits.axes_grid1.inset_locator import inset_axes
                axins2 = inset_axes(axes[1], width="100%", height="100%",                                                             # what fraction of the bounding box to take up
                                    bbox_to_anchor=(0.25, -0.2, 0.5, 0.07), bbox_transform=axes[1].transAxes)                         # x start y start x width y height    
                cb = fig.colorbar(ifg, cax = axins2, orientation = 'horizontal')
                cb.set_label("LOS displacement (m)")
                
                
            else:                                                                       # else not on a point
                pass
        else:                                                                           # else not in the axes
            pass
        fig.canvas.draw_idle()
    
    

    fig.canvas.mpl_connect("button_press_event", click)                                # connect the figure and the function, note that done when the mouse clicks.  
    axes[0].set_xlabel('Primary acquisition')
    axes[0].set_ylabel('Secondary acquisition')
    # Set the tick labels to be in dates and not pixels of the giant figure.  
    xtick_labels = []
    for tick in axes[0].get_xticks():
        try:                                                                                                    # ticks can extend past data 
            xtick_labels.append(tbaseline_info['acq_dates'][int((tick / ifg_resolution) * acq_spacing)])        
        except:
            xtick_labels.append('')
    axes[0].set_xticklabels(xtick_labels, rotation = 315, ha = 'left')
    ytick_labels = []
    for tick in axes[0].get_yticks():
        try:                                                                                                    # ticks can extend past data 
            ytick_labels.append(tbaseline_info['acq_dates'][int((tick / ifg_resolution) * acq_spacing)])        
        except:
            ytick_labels.append('')
    axes[0].set_yticklabels(ytick_labels)
    
# #%%
# import matplotlib.pyplot as plt
# import numpy as np

# def hover(event):
#     if event.inaxes == ax:                                                       # determine if the mouse is in the axes
#         print('here1')
#         if im.contains(event):                                              # cont is a boolean of if hoving on point, ind is a dictionary about the point being hovered over.  Note that two or more points can be in this.  
            
#             print(f"x:{event.xdata} y:{event.ydata}")
#         else:                                                                       # else not on a point
#             print('not contains')
#     else:                                                                           # else not in the axes
#         print("not in ax")
    
    
    
# fig, ax = plt.subplots()
# im = ax.imshow(np.random.rand(20,20))

# fig.canvas.mpl_connect("motion_notify_event", hover)                                # connect the figure and the function.  