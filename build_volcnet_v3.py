#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  3 11:26:25 2025

@author: matthew

Pseudo code

- open the xslx


- how to link labels and volcanoes?  names perhaps (or numbers)

"""

print("Started")

import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import sys
import os
import pdb
from pathlib import Path
import pandas as pd
from pathlib import PurePosixPath
from datetime import datetime
import pickle

#%%


    

def name_to_comet_frame_name(comet_frame_names, volc_name):
    """ Given a dict of regions and their frames (comet_frame_names)
    and a volcano name, find the frames for that volcano.  
    
    Inputs:
        comet_frame_names | dict of lists | e.g. south america: 
            [ 'puyuhuapi_164A_13341_131313', 'quetrupillan_083D_12832_131313', 
             'quetrupillan_083D_13027_131313',]
        volc_nae | str | e.g. puyuhuapi (camel case)
        
    Returns:
        volc_names | list of tuples | frame and region
    """
    
    volc_frames = []
        
    for region, frames in comet_frame_names.items():
        for frame in frames:
            if volc_name in frame:
                volc_frames.append((region, frame))
                
    return volc_frames

#%%

def check_or_download(jasmin_local_dir, region, file, filt = True):
    """
    
    Inputs:
        
        
    """
    
    def download_from_jasmin(
            leeds_user, leeds_server, local_path,
            jasmin_user, jasmin_server, remote_path
            ):
        """ Download a file from Jasmin to a machine that doesn't have
        access to Jasmin (hence the proxy jump)
        
        Inputs
        
        """
        
        import subprocess
        
        # make sure that all local directories exist:
        local_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Build an rsync command that uses ProxyJump for the SSH connection
        cmd = (
            f"rsync -av --progress "
            f"-e 'ssh -o ProxyJump={leeds_user}@{leeds_server}' "
            f"{jasmin_user}@{jasmin_server}:{remote_path} {local_path}"
        )
        
        result = subprocess.run(cmd, shell=True, check=True)
        
        print("STDOUT:", result.stdout)
        print("STDERR:", result.stderr)
        print("Return Code:", result.returncode)
    
    # shouldn't need to change (unless something changes on jasmin)
    from pathlib import PurePosixPath
    licsbas_path = PurePosixPath(
        "/gws/nopw/j04/nceo_geohazards_vol1/projects""/LiCS/volc-portal/"
        "processing_output/licsbas"
        )
    
    # file doens't have extension, add it and the filered identifier.  
    if filt:
        file = f"{file}_filt.json"
    else:
        file = f"{file}.json"
        
    local_path = jasmin_local_dir / region / file
    
    # Force to build a linux path as Jasmin uses this.  
    remote_path = licsbas_path / region / file

    
    # Check if the file exists
    if os.path.isfile(local_path):
        print(f"The file '{local_path}' already exists locally.")
    else:
        print(
            f"The file '{local_path}' does not exist locally so we will "
            f"attempt to download it from Jasmin"
            )

        download_from_jasmin(
            'earmgad', 'foe-linux-01.leeds.ac.uk', local_path,
            'mgaddes', 'xfer-vm-01.jasmin.ac.uk', remote_path
            )
    
    return local_path

def extract_dates(row: pd.Series):
    """
    Extracts start and end dates from a pandas Series representing a volcano eruption row.
    
    The function expects the following columns to be present:
        'Start Year', 'Start Month', 'Start Day'
        'End Year', 'End Month', 'End Day'
    
    Returns:
        A tuple (start_date, end_date) as datetime objects.
        If any conversion fails (e.g., due to missing values), that date is set to None.
    """
    from datetime import datetime
    import pandas as pd
    
    try:
        # Cast the year, month, day values to integers.
        start_year = int(float(row['Start Year']))
        start_month = int(float(row['Start Month']))
        start_day = int(float(row['Start Day']))
        start_date = datetime(start_year, start_month, start_day)
    except (ValueError, TypeError):
        start_date = None

    try:
        end_year = int(float(row['End Year']))
        end_month = int(float(row['End Month']))
        end_day = int(float(row['End Day']))
        end_date = datetime(end_year, end_month, end_day)
    except (ValueError, TypeError):
        end_date = None
        
    return start_date, end_date




def all_eruptions_one_volc(eruptions):
    """ Compiles a list of all the eruptions that overlapped with a
    Sentinel-1 time series.  
    
    Inputs:
        eruptions | data frame | each eruption for that volcano is a row
        acq_dates | list of strings | ['yyyymmdd ']
        
    Returns:
        eruption_dates | list of strings | eruptions in the form 
                                            ['yyyymmdd_yyyymmdd']
                                            
    History:
        2025_04_10 | MEG | Written
    """
    
    eruption_dates = []
    # iterate through the eruptions
    for eruption_n, eruption in eruptions.iterrows():
        # extract the erutpion date from the eruption series (dataframe row)
        erup_start, erup_end = extract_dates(eruption)
        
        # convert to a string
        eruption_dates.append(
            f"{erup_start.strftime('%Y%m%d')}_"
            f"{erup_end.strftime('%Y%m%d')}"
            )
            
    return eruption_dates


#%% currently using simpler one above that doesn't echck the intersection.  


# def all_eruptions_one_volc(eruptions, acq_dates):
#     """ Compiles a list of all the eruptions that overlapped with a
#     Sentinel-1 time series.  
    
#     Inputs:
#         eruptions | data frame | each eruption for that volcano is a row
#         acq_dates | list of strings | ['yyyymmdd ']
        
#     Returns:
#         eruption_dates | list of strings | eruptions in the form 
#                                             ['yyyymmdd_yyyymmdd']
                                            
#     History:
#         2025_04_10 | MEG | Written
#     """
    
#     return start_date, end_date



#     def check_overlap(start1: datetime, end1: datetime,
#                       start2: datetime, end2: datetime) -> (bool, int):
#         """
#         Check if two time intervals overlap.
        
#         Args:
#           start1, end1: Start and end datetimes of the first interval.
#           start2, end2: Start and end datetimes of the second interval.
        
#         Returns:
#           A tuple (overlap, gap_days) where:
#             - overlap is True if the intervals overlap (including touching edges); otherwise False.
#             - gap_days is the number of days between the intervals if they do not overlap; otherwise 0.
#         """
#         # Calculate the later of the two start dates and the earlier of the two end dates
#         latest_start = max(start1, start2)
#         earliest_end = min(end1, end2)
        
#         # If the latest start is on or before the earliest end, the intervals overlap (or touch)
#         if latest_start <= earliest_end:
#             return True, 0
#         else:
#             # Determine the gap in days between the two intervals.
#             # If the first interval ends before the second starts:
#             if end1 < start2:
#                 gap = (start2 - end1).days
#             # Or if the second interval ends before the first starts:
#             else:
#                 gap = (start1 - end2).days
            
#             return False, gap
    
#     eruption_dates = []
#     # iterate through the eruptions
#     for eruption_n, eruption in eruptions.iterrows():
#         # extract the erutpion date from the eruption series (dataframe row)
#         erup_start, erup_end = extract_dates(eruption)
        
#         # determine if hte time series and eruption overlap
#         overlap, gap_days = check_overlap(
#             datetime.strptime(acq_dates[0],'%Y%m%d'),
#             datetime.strptime(acq_dates[-1],'%Y%m%d'),
#             erup_start, erup_end
#             )
        
#         if overlap:
#             eruption_dates.append(
#                 f"{erup_start.strftime('%Y%m%d')}_"
#                 f"{erup_end.strftime('%Y%m%d')}"
#                 )
            
#     return eruption_dates


#%%


# def plot_volcano_results(image: np.ndarray, date_strings: list, 
#                          intervals: list, png_path = None, title = None):
#     """
#     Plots:
#       1. The last image from a rank-3 image array (occupying 4/5 of the vertical space)
#       2. The time series of the pixel with the largest absolute value (occupying 1/5 of the vertical space)
#       3. For each date interval (in the form 'yyyymmdd_yyyymmdd'), a horizontal line is drawn on the time
#          series plot with dots at both ends.
    
#     The legend for the time series plot is placed outside (to the right) of the axis.
    
#     Parameters:
#       image : np.ndarray
#         A rank-3 image array with shape (T, H, W) where T is time.
#       date_strings : list of str
#         A list of dates corresponding to each image frame, in the format 'yyyymmdd'.
#       intervals : list of str
#         A list of date intervals formatted as 'yyyymmdd_yyyymmdd'.
#     """
#     # Convert the date strings into datetime objects
#     time_dates = [datetime.strptime(date, "%Y%m%d") for date in date_strings]
    
#     # Create a figure with two subplots:
#     # - The top axis (ax_img) takes 4/5 of the vertical space to display the image.
#     # - The bottom axis (ax_ts) takes 1/5 of the space for the time series.
#     fig, (ax_img, ax_ts) = plt.subplots(2, 1, figsize=(12, 10), 
#                                           gridspec_kw={'height_ratios': [4, 1]})
    
#     # ----- Top Panel: Show the Last Image -----
#     last_image = image[-1]
#     im = ax_img.imshow(last_image)
#     if title is not None:
#         ax_img.set_title(title)
#     fig.colorbar(im, ax=ax_img)
    
#     # ----- Identify the Pixel of Interest -----
#     # For each pixel location (i,j), compute its maximum absolute value over time.
#     max_abs_per_pixel = np.max(np.abs(image), axis=0)
#     # Get the (row, column) index of the pixel with the largest absolute value.
#     pixel_idx = np.unravel_index(np.argmax(max_abs_per_pixel), max_abs_per_pixel.shape)
#     print("Pixel with max absolute value is at:", pixel_idx)
    
#     # Extract the time series for that pixel over all time steps.
#     pixel_series = image[:, pixel_idx[0], pixel_idx[1]]
    
#     # ----- Bottom Panel: Plot the Time Series and Intervals -----
#     ax_ts.plot(
#         time_dates, pixel_series, marker='o', linestyle='-', c = 'k',
#         label='Max. deformation'
#         )
#     ax_ts.set_xlabel("Date")
#     ax_ts.set_ylabel("Pixel Value")
#     ax_ts.set_title("Time Series of Pixel with Largest Absolute Value")
    
#     # Determine a baseline y-value for the intervals.
#     # We choose a base value a bit below the minimum of the pixel_series.
#     y_min = np.min(pixel_series)
#     y_range = np.ptp(pixel_series)  # peak-to-peak value (max - min)
#     offset = 0.1 * y_range if y_range != 0 else 1.0
#     base_y = y_min - offset

#     # Use a discrete color palette (tab10) for the intervals.
#     colors = plt.cm.tab10.colors

#     # Plot each interval as a horizontal line with dots at both ends.
#     for idx, interval in enumerate(intervals):
#         # Split the interval string "yyyymmdd_yyyymmdd" into start and end dates.
#         start_str, end_str = interval.split('_')
#         start_date = datetime.strptime(start_str, "%Y%m%d")
#         end_date = datetime.strptime(end_str, "%Y%m%d")
#         # Offset each interval vertically slightly to avoid overlap.
#         y_interval = base_y - idx * (offset * 0.5)
#         color = colors[idx % len(colors)]
#         # Draw the horizontal line representing the interval.
#         ax_ts.hlines(y=y_interval, xmin=start_date, xmax=end_date, color=color,
#                      linewidth=4, label=f"Eruption {interval}")
#         # Draw dots at both the start and end of the interval.
#         ax_ts.plot([start_date, end_date], [y_interval, y_interval], 'o', 
#                    color=color, markersize=8)
    
#     # Place the legend to the right of the time series axis.
#     ax_ts.legend(loc='center left', bbox_to_anchor=(1, 0.5))
#     # Adjust the figure to make room for the legend.
#     fig.subplots_adjust(right=0.8)
    
#     plt.tight_layout()
#     plt.show()
#     if png_path is not None:
#         fig.savefig(png_path)
#         plt.close(fig)
        

def plot_volcano_results(image: np.ndarray, dem_image: np.ndarray, date_strings: list, 
                         intervals: list, png_path=None, title=None):
    """
    Plots:
      1. The last image from a rank-3 image array (occupying the left half of the top row)
      2. The DEM image on the right half of the top row (using the 'terrain' colormap)
      3. The time series of the pixel with the largest absolute value (occupying the full width of the bottom row)
      4. For each date interval (formatted as 'yyyymmdd_yyyymmdd'), a horizontal line is drawn on the time
         series plot with dots at both ends.
    
    The legend for the time series plot is placed outside (to the right) of its axis.
    
    Parameters:
      image : np.ndarray
        A rank-3 image array with shape (T, H, W) where T is time.
      dem_image : np.ndarray
        A 2D array (H, W) representing the DEM image (assumed to have the same shape as each frame in `image`).
      date_strings : list of str
        A list of dates corresponding to each image frame, in the format 'yyyymmdd'.
      intervals : list of str
        A list of date intervals formatted as 'yyyymmdd_yyyymmdd'.
      png_path : str, optional
        If provided, the figure will be saved to this path.
      title : str, optional
        Title for the left image panel.
    """
    
    def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
        """ Take a colorbar and crop it.  Useful for removing blue parts of "terrain""
        """
        import matplotlib.colors as colors
        import numpy as np
        
        new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
        return new_cmap 
    
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec


    
    # Convert date strings into datetime objects.
    time_dates = [datetime.strptime(date, "%Y%m%d") for date in date_strings]
    
    # Create the figure and GridSpec layout:
    # - 2 rows and 2 columns, with the top row (images) having a height ratio of 4 and bottom row (time series) 1.
    # - The time series spans both columns.
    fig = plt.figure(figsize=(12, 10))
    gs = GridSpec(nrows=2, ncols=2, height_ratios=[4, 1], width_ratios=[1, 1], wspace=0.3)
    
    # Top row: left panel for the primary image and right panel for the DEM.
    ax_img = fig.add_subplot(gs[0, 0])
    ax_dem = fig.add_subplot(gs[0, 1])
    # Bottom row: time series across both columns.
    ax_ts = fig.add_subplot(gs[1, :])
    
    # ----- Top Left Panel: Plot the Last Image -----
    # Take the last frame from the rank-3 image array.
    last_image = image[-1]
    im = ax_img.imshow(last_image)
    if title is not None:
        ax_img.set_title(title)
    cbar = fig.colorbar(im, ax=ax_img)
    cbar.set_label("LOS Displacement (m)")
    
    # ----- Top Right Panel: Plot the DEM Image -----
    terrain_cmap = truncate_colormap(
        plt.get_cmap('terrain') , 0.2, 1
        )
    dem_im = ax_dem.imshow(dem_image, cmap=terrain_cmap)
    ax_dem.set_title("DEM")
    cbar_dem = fig.colorbar(dem_im, ax=ax_dem)
    cbar_dem.set_label("Elevation (m)")
    
    # ----- Identify the Pixel of Interest (for the time series) -----
    # For each pixel location (i,j), compute its maximum absolute value over time.
    max_abs_per_pixel = np.max(np.abs(image), axis=0)
    # Get the (row, column) index of the pixel with the largest absolute value.
    pixel_idx = np.unravel_index(np.argmax(max_abs_per_pixel), max_abs_per_pixel.shape)
    print("Pixel with max absolute value is at:", pixel_idx)
    # Extract the time series for that pixel over all time steps.
    pixel_series = image[:, pixel_idx[0], pixel_idx[1]]
    
    # ----- Bottom Panel: Plot the Time Series and Intervals -----
    ax_ts.plot(
        time_dates, pixel_series, marker='o', linestyle='-', c='k',
        label='Max. deformation'
    )
    ax_ts.set_xlabel("Date")
    ax_ts.set_ylabel("LOS Displacement (m)")
    ax_ts.set_title("Time Series of Pixel with Largest Absolute Value")
    
    # Determine a baseline y-value for the intervals.
    y_min = np.min(pixel_series)
    y_range = np.ptp(pixel_series)  # peak-to-peak value (max - min)
    offset = 0.1 * y_range if y_range != 0 else 1.0
    base_y = y_min - offset

    # Use a discrete color palette (tab10) for the intervals.
    colors = plt.cm.tab10.colors
    # Plot each interval as a horizontal line with dots at both ends.
    for idx, interval in enumerate(intervals):
        # Split the interval string "yyyymmdd_yyyymmdd" into start and end dates.
        start_str, end_str = interval.split('_')
        start_date = datetime.strptime(start_str, "%Y%m%d")
        end_date = datetime.strptime(end_str, "%Y%m%d")
        # Offset the interval vertically slightly to avoid overlap.
        y_interval = base_y - idx * (offset * 0.5)
        color = colors[idx % len(colors)]
        # Draw the horizontal line representing the interval.
        ax_ts.hlines(y=y_interval, xmin=start_date, xmax=end_date, color=color,
                     linewidth=4, label=f"Eruption {interval}")
        # Draw dots at both the start and end of the interval.
        ax_ts.plot([start_date, end_date], [y_interval, y_interval], 'o', 
                   color=color, markersize=8)
    
    # Place the legend to the right of the time series axis.
    ax_ts.legend(loc='center left', bbox_to_anchor=(1.02, 0.5))
    
    # Adjust the layout to leave space for the external legend without shifting the top panels.
    plt.tight_layout(rect=[0, 0, 0.95, 1])
    plt.show()
    
    if png_path is not None:
        fig.savefig(png_path)
        plt.close(fig)
            

#%%


licsalert_pkg_dir  = Path(
    "/home/matthew/university_work/03_automatic_detection_algorithm/"
    "06_LiCSAlert/00_LiCSAlert_GitHub"
    )

# add to path if not already
if str(licsalert_pkg_dir) not in sys.path:                              
    sys.path.append(str(licsalert_pkg_dir))                            

import licsalert
from licsalert.jasmin_tools import open_comet_frame_files
from licsalert.data_importing import LiCSBAS_json_to_LiCSAlert


#%% inputs

# reord of eruptions
smithsonian_xslx_path = Path("GVP_Eruption_Search_Result_202504030559.xlsx")

# volc id to jasmin_names
mappings_path = Path("mappings.xlsx")

# all comet volcano frames by region
volc_frame_names = Path("./comet_volcano_frames")

# licsbas path (on Jasmin_)
jasmin_licsbas_dir = PurePosixPath(
    "/gws/nopw/j04/nceo_geohazards_vol1/projects/LiCS/" 
    "volc-portal/processing_output/licsbas"
    )

# local copy of Jasmin files
jasmin_local_dir = Path('./jasmin_clones')


# labelled data and png samples stored here
volcnet_outdir = Path('./volcnet_labelled_data')

# delete local .json files after use (to save space)
cleanup_json = True

# overlap_required


# supress figures
plt.switch_backend('Agg')

#%% Open Smithsonian data, mappings (volc_id to name), 
# and comet frame names

eruptions_smithsonian = pd.read_excel(smithsonian_xslx_path, skiprows=1)

# get column names
# eruptions_smithsonian.columns

#get only the eruptions since S1 started.  
eruptions_s1 = eruptions_smithsonian[
    eruptions_smithsonian['Start Year'] > 2015
    ]
    

# also get the mappings from volcid to name
mappings = pd.read_excel(mappings_path)


# comet frame names by region 
comet_frame_names = open_comet_frame_files(volc_frame_names)


#%%  Need to collect all eruptions by volcano type.  

# Group the DataFrame by "Volcano Name", and sort alphabetically
eruptions_by_volc = sorted(
    eruptions_s1.groupby('Volcano Name'), key=lambda x: x[0]
    )


       
                
for volcano_name, eruptions in eruptions_by_volc:
    print(f"Volcano Name: {volcano_name}")

    # Get eruption dates for this volcano.
    try:
        eruption_dates = all_eruptions_one_volc(eruptions)
    except Exception as e:
        print(f"Error getting eruption dates for {volcano_name}: {e}")
        continue

    # Get the volcano number (there should be exactly one).
    try:
        unique_numbers = eruptions['Volcano Number'].unique()
        if len(unique_numbers) != 1:
            raise Exception(
                "Multiple volcano numbers found for a single volcano."
                )
        volc_n = unique_numbers[0]
    except Exception as e:
        print(f"Error with Volcano Number for {volcano_name}: {e}")
        continue

    # Convert volcano number to a jasmin name using the mappings.
    try:
        mapping = mappings[mappings.iloc[:, 1] == volc_n]['jasmin_name']
        if mapping.empty:
            raise Exception("Volcano number not found in COMET mappings.")
        jasmin_name = mapping.values[0]
    except Exception as e:
        print(
            f"Error converting Volcano Number to jasmin name for"
            f" {volcano_name}: {e}"
            )
        continue

    # Get the comet frames for that volcano.
    try:
        frames = name_to_comet_frame_name(comet_frame_names, jasmin_name)
    except Exception as e:
        print(f"Error retrieving COMET frames for {volcano_name}: {e}")
        continue

    # Process each frame.
    for frame in frames:
        try:
            json_path = check_or_download(
                jasmin_local_dir, region=frame[0], file=frame[1]
                )
        except Exception as e:
            print(
                f"Error in check_or_download for {volcano_name}, frame "
                "{frame}: {e}"
                )
            continue

        try:
            outputs = LiCSBAS_json_to_LiCSAlert(
                json_path, crop_side_length=None, mask_type='nan_variable'
            )
            _, displacement_r3, tbaseline_info, ref_xy, json_time = outputs

            # Create DEM water mask.
            displacement_r3['water_mask'] = (displacement_r3['dem'] < 1.e-19)
        except Exception as e:
            print(
                f"Error opening the json file or creating the water mask "
                "for {volcano_name}, frame {frame}: {e}"
                )
            continue

        try:
            # Ensure the output directory exists.
            (volcnet_outdir / frame[0]).mkdir(parents=True, exist_ok=True)
        except Exception as e:
            print(f"Error creating output directory for {volcano_name}, frame {frame}: {e}")
            continue

        try:
            # Create the PNG plot.
            plot_volcano_results(
                displacement_r3['cumulative'],
                ma.array(
                    displacement_r3['dem'], 
                    mask = displacement_r3['water_mask']
                    ),
                tbaseline_info['acq_dates'],
                eruption_dates, 
                png_path=volcnet_outdir / frame[0] / f"{json_path.stem}.png",
                title=json_path.stem
            )
            
        except Exception as e:
            print(
                f"Error plotting volcano results for {volcano_name}, frame "
                "{frame}: {e}"
                )
            continue

        try:
            # Remove unneeded keys.
            del displacement_r3['cumulative']
            del tbaseline_info['ifg_dates']
            del tbaseline_info['baselines']
            del tbaseline_info['baselines_cumulative']
            
            # Save the processed data to a pickle file.
            pkl_name = volcnet_outdir / frame[0] / f"{json_path.stem}.pkl"
            with open(pkl_name, 'wb') as f:
                pickle.dump(displacement_r3, f)
                pickle.dump(tbaseline_info, f)
                pickle.dump(eruption_dates, f)
        except Exception as e:
            print(
                f"Error saving pickle for {volcano_name}, frame {frame}: {e}"
                )
            continue

        if cleanup_json:
            try:
                print(f"Deleting {json_path} to save disk space")
                json_path.unlink()
            except Exception as e:
                print(
                    f"Error deleting JSON file {json_path} for {volcano_name} "
                    ", frame {frame}: {e}"
                    )
                continue
        
            

#%% Filter after making all of them




#%%

sys.exit()

