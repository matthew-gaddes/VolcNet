#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 19 12:59:21 2025

@author: earcno
"""

import csv
import numpy as np
from collections import Counter
import matplotlib.pyplot as plt
from datetime import datetime
import matplotlib.dates as mdates


# Read the .csv file
file_path = './GVP_Eruption_Search_Result_202503190715(Eruption List).csv'
data = []

with open(file_path, newline='', encoding='latin1') as csvfile:
    reader = csv.reader(csvfile)
    next(reader)  # Skip the first line (header)
    next(reader)  # Skip the second line
    for row in reader:
        data.append(row)

data = np.array(data)

# Eruptions from 2015: 336
eruptions_from_2015 = data[data[:, 8].astype(int) >= 2015]

volcano_index=eruptions_from_2015[:,0]
names = eruptions_from_2015[:, 1]

# Count unique names and their repetitions in column 2 (index 1)
name_counts = Counter(names)

# Create a bar graph for each name and count
plt.figure(figsize=(10, 6))
plt.bar(name_counts.keys(), name_counts.values())
plt.xlabel('Names')
plt.ylabel('Counts')
plt.title('Counts of Unique Names in Filtered Data')
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig('name_counts_bar_graph.png')
plt.show()

# Count the number of unique names
number_of_unique_names = len(name_counts)
print(f"Number of unique names: {number_of_unique_names}")

# Plot the location of each volcano on a map using latitude and longitude in columns 23 and 24 (indexes 22 and 23)
lat_column_index = 22
lon_column_index = 23
latitudes = eruptions_from_2015 [:, lat_column_index].astype(float)
longitudes = eruptions_from_2015 [:, lon_column_index].astype(float)

# Assign colors based on counts
unique_names = list(name_counts.keys())
colors = plt.cm.rainbow(np.linspace(0, 1, len(unique_names)))
name_to_color = {name: colors[i] for i, name in enumerate(unique_names)}
point_colors = [name_to_color[name] for name in names]


# Assign colors based on counts
unique_names = list(name_counts.keys())
counts = list(name_counts.values())
colors = plt.cm.rainbow(np.linspace(0, 1, max(counts)))
count_to_color = {count: colors[count-1] for count in counts}
point_colors = [count_to_color[name_counts[name]] for name in names]

# Create a plot with the world map and volcano locations
fig, ax = plt.subplots(figsize=(15, 10))
ax.scatter(longitudes, latitudes, marker='o', c=point_colors, s=50)

# Create a legend for the colors based on counts
unique_counts = sorted(set(counts))
handles = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=count_to_color[count], markersize=10) for count in unique_counts]
ax.legend(handles, unique_counts, title="Counts", bbox_to_anchor=(1.05, 1), loc='upper left')

ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
ax.set_title('Locations of Volcanoes')
plt.grid(True)
plt.tight_layout()
plt.savefig('volcano_locations_map_with_colors_and_legend.png')
plt.show()

# Extract the first and second dates from the specified columns
first_dates = eruptions_from_2015[:, [8, 10, 12]]
second_dates = eruptions_from_2015[:, [16, 18, 20]]

# Calculate the number of days between the first date and the second date for each row
days_between = []
for i in range(len(first_dates)):
    first_date_str = f"{int(first_dates[i, 0])}-{int(first_dates[i, 1]):02d}-{int(first_dates[i, 2]):02d}"
    second_date_str = f"{int(second_dates[i, 0])}-{int(second_dates[i, 1]):02d}-{int(second_dates[i, 2]):02d}"
    
    first_date = datetime.strptime(first_date_str, "%Y-%m-%d")
    second_date = datetime.strptime(second_date_str, "%Y-%m-%d")
    
    delta = (second_date - first_date).days
    days_between.append(delta)

# Display the number of days between the first date and the second date for each row
for i, days in enumerate(days_between):
    print(f"Row {i+1}: {days} days between the first date and the second date")
    
# Convert days to months (approximate)
months_between = [days / 30 for days in days_between]

# Create a bar plot of the 'months_between' data
fig, ax = plt.subplots(figsize=(10, 6))
bars = ax.bar(range(len(months_between)), months_between, color='b')
plt.xlabel('Eruption')
plt.ylabel('Duration of eruption [months]')
plt.title('Duration of eruptions from 2015')
plt.grid(True)

# Add interactivity to display name and dates on click
def on_click(event):
    for bar in bars:
        if bar.contains(event)[0]:
            index = bars.index(bar)
            name = data[index, 1]
            first_date_str = f"{int(first_dates[index, 0])}-{int(first_dates[index, 1]):02d}-{int(first_dates[index, 2]):02d}"
            second_date_str = f"{int(second_dates[index, 0])}-{int(second_dates[index, 1]):02d}-{int(second_dates[index, 2]):02d}"
            print(f"Name: {name}, First Date: {first_date_str}, Second Date: {second_date_str}")

fig.canvas.mpl_connect('button_press_event', on_click)

plt.show()

## Figure to show the duration of each eruption for each volcano

# Convert dates to datetime objects
first_dates_dt = [datetime(int(date[0]), int(date[1]), int(date[2])) for date in first_dates if date[0] and date[1] and date[2]]
second_dates_dt = [datetime(int(date[0]), int(date[1]), int(date[2])) for date in second_dates if date[0] and date[1] and date[2]]

# Extract volcano names from column 2
volcano_names = eruptions_from_2015[:, 1]

# Create a plot with the y-axis containing the name of the volcano and the x-axis containing the dates
fig, ax = plt.subplots(figsize=(10, 6))

# Plot bars between the first date and the second date for each volcano
for i, (start_date, end_date) in enumerate(zip(first_dates_dt, second_dates_dt)):
    ax.barh(volcano_names[i], (end_date - start_date).days, left=start_date)

plt.xlabel('Dates')
plt.ylabel('Volcano Names')
plt.title('Volcano Eruption Duration')
plt.grid(True)
plt.tight_layout()
plt.savefig('volcano_eruption_duration.png')
plt.show()

# Divide the data into 3 different plots

# Convert dates to datetime objects, handling missing values
first_dates_dt = []
second_dates_dt = []
for date in first_dates:
    if date[0] and date[1] and date[2]:
        first_dates_dt.append(datetime(int(date[0]), int(date[1]), int(date[2])))
    else:
        first_dates_dt.append(None)

for date in second_dates:
    if date[0] and date[1] and date[2]:
        second_dates_dt.append(datetime(int(date[0]), int(date[1]), int(date[2])))
    else:
        second_dates_dt.append(None)



# Divide the data into 5 different plots
num_plots = 5
volcanoes_per_plot = len(volcano_names) // num_plots

for i in range(num_plots):
    fig, ax = plt.subplots(figsize=(10, 6))
    
    start_idx = i * volcanoes_per_plot
    end_idx = (i + 1) * volcanoes_per_plot if i < num_plots - 1 else len(volcano_names)
    
    # Plot bars between the first date and the second date for each volcano in the current plot
    for j in range(start_idx, end_idx):
        ax.barh(volcano_names[j], (second_dates_dt[j] - first_dates_dt[j]).days, left=first_dates_dt[j])
    
    plt.xlabel('Dates')
    plt.ylabel('Volcano Names')
    plt.title(f'Volcano Eruption Duration (Plot {i + 1})')
    plt.grid(True)
    
    # Format x-axis labels as year-month-day
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
    
    plt.tight_layout()
    plt.savefig(f'volcano_eruption_duration_plot_{i + 1}.png')
    plt.show()