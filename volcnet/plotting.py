#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  3 12:33:40 2022

@author: matthew
"""

import pdb

def volcnet_ts_visualiser(displacement_r3, tbaseline_info, persistent_defs, transient_defs, 
                          acq_spacing = 5, ifg_resolution = 20, figsize_height = 10, labelling_function = None):
    """Visualise all the possible interferograms that can be formed from a VolcNet time series.  
    
    Inputs:
        displacement_r3 | dict | rank 3 data (i.e. time x y x x), contains:
                                    cumulative | rank 3 ma | cumulative displacement at each acquisiiton date, relative to the first acquisition.  
                                    mask | rank 2 | 1 where masked due to water or incoherence.  
                                    dem | rank 2 ma | DEM
                                    lons | rank 2 | longitudes for each pixel.  
                                    lats | rank 2 | latitudes for each pixel.  
        tbaselines_info | dict | contains:
                                    acq_dates | list | acquisitions, in form YYYYMMDD
                                    baselines_cumulative | rank 1 array | days since first acquistiion for each acquisition.  
                                    
        persisent_defs | list of dicts | Persistent deformations (e.g inflation of years of months).  E.g.:
                                                                        {'def_lon_west': -91.17,
                                                                         'def_lon_east': -91.11,
                                                                         'def_lat_south': -0.84,
                                                                         'def_lat_north': -0.8,
                                                                         'def_episode_start': 20141213,
                                                                         'def_episode_stop': 20170706,
                                                                         'source': 'sill',
                                                                         'def_rate': 0.4}
                                                                        
        transient_defs | list of dicts | Transient deformations. e.g. an eruption that is imaged between two acquisitions.   E.g.:
                                                                            {'def_lon_west': -91.28,
                                                                             'def_lon_east': -91.11,
                                                                             'def_lat_south': -0.84,
                                                                             'def_lat_north': 0.7,
                                                                             'def_episode_start': 20180526,
                                                                             'def_episode_stop': 20180607,
                                                                             'source': 'sill',
                                                                             'def_magnitude': 0.7}

                                                                        
        acq_spacing | int | every acq_spacing interferogram is shown.  Larger number makes a smaller (and more manageable) figure.  
        ifg_resolution | int | all the interferograms (spaced every acq_spacing) will be displayed as ifg_resolution x ifg_resolution thumbnails.  Note that aspect is destroyed by this, as thumbnails must be square (to ensure x and y time axis are equal)
        figsize_height | int | figure height in inches.  Width is double this.  
        
    Returns:
        Figure
        
    History:
        2022_05_03 | MEG | Written.  
    """

    import numpy as np
    import numpy.ma as ma
    from datetime import datetime
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    import matplotlib.gridspec as gridspec
    import matplotlib.ticker as mticker
    
    from volcnet.aux import ll_2_pixel
    
    def click(event):
        if event.inaxes == ax_all_ifgs:                                                                    # determine if the mouse is in the axes on the left
            if all_ifgs.contains(event):                                                                   # cont is a boolean of if hoving on point, ind is a dictionary about the point being hovered over.  Note that two or more points can be in this.  
                
                # 0: clear axes from the previous time we drew a single ifg.  
                try:
                    for ax in fig.axes[4:]:                                         
                        ax.remove()                                                                        # possibly remove the colobar axis so it can be drawn again     
                except:
                    pass
                ax_1_ifg.clear()                                                                            # and clear the axes to get rid of imshow and plot (for the location information)
                
                # 1: Make the ifg that the mouse is hovering on                    
                acq_primary = int(event.xdata / ifg_resolution) * acq_spacing                                                               # xdata is in pixels of the big image, figure, convert to which acquisition number this is.   
                acq_secondary = int(event.ydata / ifg_resolution) * acq_spacing                                                             # and in y 
                #print(f"x:{acq_primary} y:{acq_secondary}")
                ifg_array = displacement_r3['cumulative'][acq_primary] - displacement_r3['cumulative'][acq_secondary]                       # make the ifg we are hovering on.  
                cmap_mid = 1 - ma.max(ifg_array)/(ma.max(ifg_array) + abs(ma.min(ifg_array)))                                               # get the ratio of the data that 0 lies at (eg if data is -15 to 5, ratio is 0.75)
                colours_cent = remappedColorMap(plt.get_cmap('coolwarm'), start=0.0, midpoint=cmap_mid, stop=1, name='ic_colours_cent')     # make the colours for plotting the ICs
                ifg = ax_1_ifg.imshow(ifg_array, cmap = colours_cent)                                                                       # draw the ifg we are hovering on.  
                
                axins2 = inset_axes(ax_1_ifg, width="100%", height="100%", loc = 'upper left',                                              # what fraction of the bounding box to take up
                                    bbox_to_anchor=(0.05, 0.95, 0.3, 0.05), bbox_transform=ax_1_ifg.transAxes)                              # x start y start x width y height    
                cb = fig.colorbar(ifg, cax = axins2, orientation = 'horizontal')
                cb.set_label("LOS displacement (m)")
                
                # 2: Add the temporal baseline to the figure.  
                primary = datetime.strptime(tbaseline_info['acq_dates'][acq_primary], '%Y%m%d')                                                     # to label with the temporal baseline, need to get acquisitions.  
                secondary = datetime.strptime(tbaseline_info['acq_dates'][acq_secondary], '%Y%m%d')
                tbaseline = (primary - secondary).days
                ifg_title = f"{tbaseline_info['acq_dates'][acq_secondary]}_{tbaseline_info['acq_dates'][acq_primary]} ({tbaseline} days)"
                ax_1_ifg.set_title(ifg_title)

                
                # 3: update the tick labels to be lon lat
                xtick_labels = []
                for tick in ax_1_ifg.get_xticks():
                    try:                                                                                                    # ticks can extend past data 
                        xtick_labels.append(np.round(displacement_r3['lons'][-1, int(tick)], 2))        
                    except:
                        xtick_labels.append('')
                ax_1_ifg.xaxis.set_major_locator(mticker.FixedLocator(ax_1_ifg.get_xticks()))
                ax_1_ifg.xaxis.set_major_formatter(mticker.FixedFormatter(xtick_labels))
                plt.setp(ax_1_ifg.get_xticklabels(), rotation=35, horizontalalignment='left')
                
                ytick_labels = []
                for tick in ax_1_ifg.get_yticks():
                    try:                                                                                                    # ticks can extend past data 
                        ytick_labels.append(np.round(displacement_r3['lats'][int(tick), 0], 2))        
                    except:
                        ytick_labels.append('')
                ax_1_ifg.yaxis.set_major_locator(mticker.FixedLocator(ax_1_ifg.get_yticks()))
                ax_1_ifg.yaxis.set_major_formatter(mticker.FixedFormatter(ytick_labels))
                
                if labelling_function is not None:
                    ifg_name = f"{tbaseline_info['acq_dates'][acq_secondary]}_{tbaseline_info['acq_dates'][acq_primary]}"            # get the name of the ifg in yyyymmdd_yyyymmdd
                    def_predicted, sources, def_location = labelling_function(ifg_name, persistent_defs, transient_defs)             # determine label.
                    ifg_title = ifg_title + f"\n Source(s): {sources} Magnitude: {def_predicted:.2f}m"                               # add label info to title, deformation is to 2dp.  
                    ax_1_ifg.set_title(ifg_title)                                                                                    # make label
                    xys_array = ll_2_pixel(def_location, displacement_r3['lons'], displacement_r3['lats'])                           # convert lon and lat of deformation to pixel number
                    ax_1_ifg.plot(xys_array[:,0], xys_array[:,1] )                                                                   # plot pixel numbers       
                    

            else:                                                                                                           # else not on a point
                pass
        else:                                                                                                               # else not in the axes
            pass
        fig.canvas.draw_idle()
    
    n_acq, ny, nx = displacement_r3['cumulative'].shape

    
    # 1: Create the large figure showing ifgs for each acquisition pair).  
    cumulative_r3_small = displacement_r3['cumulative'][::acq_spacing,
                                                        np.array((np.linspace(0, ny-1, ifg_resolution)), dtype = int), :]            # downsample in time and y
    cumulative_r3_small = cumulative_r3_small[:,:, np.array((np.linspace(0, nx-1, ifg_resolution)), dtype = int)]                    # downsample in x
    n_acq_plot, ny_plot, nx_plot = cumulative_r3_small.shape                                                                         # 
    figure = ma.zeros((n_acq_plot * ifg_resolution, n_acq_plot * ifg_resolution))                                                    # create a giant ma array of zeros to hold all the one channel images.  

    for row_n in range(n_acq_plot):                                                                                                  # loop down one row.  
        print(f"row {row_n}")
        cumulative_r3_small -= cumulative_r3_small[row_n,]                                                                           # make the time series relative to this acquisition number.  
        for col_n in range(n_acq_plot):
            if col_n == row_n:                                                                                                       # if on the diagonal, make a square of zeros so it stands out
                mini_ifg = ma.array(np.zeros((ny_plot, nx_plot)), mask = np.zeros((ny_plot, nx_plot)))
            else:
                mini_ifg = cumulative_r3_small[col_n,]                                                                              # for all others just get the ifg.  
            figure[row_n*ifg_resolution : (row_n + 1) *ifg_resolution,  
                   col_n*ifg_resolution : (col_n + 1) *ifg_resolution,] = mini_ifg                                                  # put the ifg in the correct place on the giant figure.              
            
    

    #  2: start the figure:
    fig = plt.figure(figsize=(2*figsize_height, figsize_height))
    grid = gridspec.GridSpec(12, 24, wspace=0.1, hspace=0.1)                                                    # 
    ax_all_ifgs = plt.Subplot(fig, grid[:12, :12])                                                                                     # left hand axes that is square.  
    ax_1_ifg = plt.Subplot(fig, grid[:8, 12:])                                                                                         # axes on right to show 1 ifg.  
    ax_labels = plt.Subplot(fig, grid[9:, 13:-1])                                                                                      # lower right axes to show labels.  
    
    # cmap_mid = 1 - ma.max(figure)/(ma.max(figure) + abs(ma.min(figure)))                                                             # get the ratio of the data that 0 lies at (eg if data is -15 to 5, ratio is 0.75)
    # figure_cmap = remappedColorMap(plt.get_cmap('coolwarm'), start=0.0, midpoint=cmap_mid, stop=1, name='ic_colours_cent')           # make the colours for plotting the ICs
    all_ifgs = ax_all_ifgs.imshow(figure)                                                                                              # plot the giant overview of tiny interferograms.  

    fig.canvas.mpl_connect("button_press_event", click)                                                                                 # connect the figure and the function, note that done when the mouse clicks.  
    ax_all_ifgs.set_xlabel('Primary acquisition')
    ax_all_ifgs.set_ylabel('Secondary acquisition')
    # Set the tick labels to be in dates and not pixels of the giant figure.  
    xtick_labels = []
    for tick in ax_all_ifgs.get_xticks():
        try:                                                                                                    # ticks can extend past data 
            xtick_labels.append(tbaseline_info['acq_dates'][int((tick / ifg_resolution) * acq_spacing)])        
        except:
            xtick_labels.append('')

    ax_all_ifgs.xaxis.set_major_locator(mticker.FixedLocator(ax_all_ifgs.get_xticks()))
    ax_all_ifgs.xaxis.set_major_formatter(mticker.FixedFormatter(xtick_labels))
    ax_all_ifgs.xaxis.tick_top()
    ax_all_ifgs.xaxis.set_label_position('top')
    plt.setp(ax_all_ifgs.get_xticklabels(), rotation=35, horizontalalignment='left')
    
    ytick_labels = []
    for tick in ax_all_ifgs.get_yticks():
        try:                                                                                                    # ticks can extend past data 
            ytick_labels.append(tbaseline_info['acq_dates'][int((tick / ifg_resolution) * acq_spacing)])        
        except:
            ytick_labels.append('')
    ax_all_ifgs.yaxis.set_major_locator(mticker.FixedLocator(ax_all_ifgs.get_yticks()))
    ax_all_ifgs.yaxis.set_major_formatter(mticker.FixedFormatter(ytick_labels))

    fig.add_subplot(ax_all_ifgs)                                                                   # add to figure
    fig.add_subplot(ax_1_ifg)                                                                   # add to figure
    
    
    # 3: Add the label (as in labelled data) information at the bottom.  
    fig.add_subplot(ax_labels)                                                                   # add to figure
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
                


########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
#################                                                                                                      #################
#################                        Reproduce from insar_tools                                                    #################
#################                                                                                                      #################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################



#%%



def xticks_every_nmonths(ax_to_update, day0_date, time_values, include_tick_labels, 
                         major_ticks_n_months = 3, minor_ticks_n_months = 1):
    """Given an axes, update the xticks so the major ones are the 1st of every n months (e.g. if every 3, would be: jan/april/july/october).  
    
    Inputs:
        ax_to_update | matplotlib axes | the axes to update.  
        day0_date | string | in form yyyymmdd
        time_values | rank 1 array | cumulative temporal baselines, e.g. np.array([6,18, 30, 36, 48])
        include_tick_labels | boolean | if True, tick labels are added to the ticks.  
        n_months | int | x ticks are very n months.  e.g. 2, 3,4,6,12 (yearly)  Funny spacings (e.g. every 5) not tested.  
    Returns:
        updates axes
    History:
        2021_09_27 | MEG | Written
        2022_02_17 | MEG | modify so can be monhtly spacing other than every 3 months.  
    """
    import pdb
    import numpy as np
    import datetime as dt
    import copy
    
    from dateutil.relativedelta import relativedelta                                                    # add 3 months and check not after end
    import matplotlib.pyplot as plt

    
    def create_tick_labels_every_nmonths(day0_date_dt, dayend_date_dt, n_months = 1):
        """ Given a spacing of every n_months, get the dates and days since the first date for ticks every n_months.  
        e.g. every month, every 6 months.  
        
        Inputs:
            day0_date_dt | datetime | date of x = 0 on axis.  
            dayend_date_dt | datetime | date of last x value. 
            n_months | int | frequency of ticks. 
            
        Returns:
            ticks | dict | contains datetimes : datetimes for each tick
                                    yyyymmdd : strings to use as labels in form yyyy/mm/dd
                                    n_day     : day number of tick.  
        History:
            2022_03_29 | MEG | Written
        """
    
        # 1: find first tick date (the first of the jan/ april/jul /oct)                        
        date_tick0 = copy.deepcopy(day0_date_dt)                                                                 # version that can be modified as we iterate through.  
        while not ( (date_tick0.day) == 1 and (date_tick0.month in (np.arange(0, 12, n_months) + 1))):           # i.e. whilst it's not the 1st of the month, and not jan/apr/jul/oct....
            date_tick0 +=  dt.timedelta(1)                                                                       # then add one day and keep going.  
    
        # 2: get all the other first of the quarters as datetimes (ie keep adding n months until we're gone past the day end date)
        ticks = {'datetimes' : [date_tick0],
                 'yyyymmdd'   : [],
                 'n_day'     : []}
       
        while ticks['datetimes'][-1] < (dayend_date_dt - relativedelta(months=+n_months)):                         # while we haven't gone past the last date (and subtract 3 months to make sure we don't go one 3 month jump too far. )
            ticks['datetimes'].append(ticks['datetimes'][-1] + relativedelta(months=+n_months))                    # append the next date which is n_months more.  
        
        # 3: work out what day number each first of the quarter is.  
        for tick_dt in ticks['datetimes']:                                                                      # loop along the list of datetimes (which are each tick) 
            ticks['yyyymmdd'].append(dt.datetime.strftime(tick_dt, "%Y/%m/%d"))                                 # as a string that can be used for the tick label (hence why include / to make more readable)
            ticks['n_day'].append((tick_dt - day0_date_dt).days)
            
        return ticks

    
    xtick_label_angle = 315                                                                              # this angle will read from top left to bottom right (i.e. at a diagonal)
    
    tick_labels_days = ax_to_update.get_xticks().tolist()                                                # get the current tick labels
    day0_date_dt = dt.datetime.strptime(day0_date, "%Y%m%d")                                             # convert the day0 date (date of day number 0) to a datetime.  
    dayend_date_dt = day0_date_dt +  dt.timedelta(int(time_values[-1]))                                  # the last time value is the number of days we have, so add this to day0 to get the end.  

    ticks_major = create_tick_labels_every_nmonths(day0_date_dt, dayend_date_dt, n_months = major_ticks_n_months)
    ticks_minor = create_tick_labels_every_nmonths(day0_date_dt, dayend_date_dt, n_months = minor_ticks_n_months)                        # these are used as the minor ticks every month.  
    

        
    # 4: Update the figure.  
    ax_to_update.set_xticks(ticks_major['n_day'])                                                                   # apply major tick labels to the figure
    ax_to_update.set_xticks(ticks_minor['n_day'], minor = True)                                                                   # apply major tick labels to the figure

    if include_tick_labels:
        ax_to_update.set_xticklabels(ticks_major['yyyymmdd'], rotation = xtick_label_angle, ha = 'left')            # update tick labels, and rotate
        plt.subplots_adjust(bottom=0.15)
        ax_to_update.set_xlabel('Date')
    else:
        ax_to_update.set_xticklabels([])                                                                    # remove any tick lables if they aren't to be used.  
    
    # add vertical lines every year.  
    for major_tick_n, datetime_majortick in enumerate(ticks_major['datetimes']):
        if datetime_majortick.month == 1:                                                                       # if it's the january tick (i.e the 1st of the year)
            ax_to_update.axvline(x = ticks_major['n_day'][major_tick_n], color='k', alpha=0.1, linestyle='--')                          
               


    
    
#%%


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