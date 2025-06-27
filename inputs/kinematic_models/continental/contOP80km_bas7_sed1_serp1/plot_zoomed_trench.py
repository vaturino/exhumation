#!/usr/bin/python3

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import sys
from matplotlib.colors import ListedColormap, BoundaryNorm

file_name = str(sys.argv[1])
cmin = float(sys.argv[2])  # min colorbar value (not needed but kept for flexibility)
cmax = float(sys.argv[3])  # max colorbar value

# Define plot limits
zmin, zmax = 690.e3, 900.e3
xmin, xmax = 2300.e3, 2750.e3

# Load file
file_path = f'text_files/{file_name}.txt'
plot_name = f'test_plots/{file_name}_zoomed.png'
file = pd.read_csv(file_path, sep=r"\s+", skiprows=2, header=None)

file.columns = ['x', 'y', 'oc', 'op', 'sed', 'opc', 'core', 'wz', 'ecl', 'serp']

# Filter data within the zoomed area
file_new = file[(file['x'] > xmin) & (file['x'] < xmax) & (file['y'] > zmin) & (file['y'] < zmax)]
file_new["x"] = file_new["x"]/1e3
file_new["y"] = (zmax - file_new["y"])/1e3
print("Read and filtered file...")

# Define terrain colors
terrain_palette = {
    "core": sns.color_palette("Paired")[1],
    "oc": "olivedrab",
    "sed": sns.color_palette("Paired")[2],
    "opc": sns.color_palette("Paired")[3],
    "op": sns.color_palette("Paired")[8],
    "serp": "darkslategray",
    "ecl": "sienna",  
    "wz": "darkorange"
}

fontsize = 16

# Convert terrain colors to a colormap
terrain_types = list(terrain_palette.keys())
terrain_colors = list(terrain_palette.values())

# Add white for background
terrain_colors.insert(0, "white")  # Background (index 0)
terrain_types.insert(0, "background")  # Corresponding label

terrain_cmap = ListedColormap(terrain_colors)
norm = BoundaryNorm(range(len(terrain_types) + 1), terrain_cmap.N)

# Identify terrain type for each point
terrain_matrix = file_new[terrain_types[1:]].values  # Exclude "background" column
terrain_indices = np.argmax(terrain_matrix, axis=1) + 1  # Shift indices by +1

# Set background (where all values are 0)
terrain_indices[np.all(terrain_matrix == 0, axis=1)] = 0  # Assign index 0 for white

# Create figure
fig, ax = plt.subplots(figsize=(10, 8))
ax.set_facecolor('white')

# Scatter plot with assigned colors
sc = ax.scatter(file_new['x'], file_new['y'], c=terrain_indices, cmap=terrain_cmap, norm=norm, s=1)

# Set plot limits
ax.set_xlim(xmin/1e3, xmax/1e3)
ax.set_ylim((zmax-zmin)/1e3, 0)
#set ticks every 100 km and the font is 12
plt.xticks(np.arange(xmin/1e3, xmax/1e3, 100), fontsize=fontsize)
plt.yticks(np.arange(0, (zmax-zmin)/1e3, 50), fontsize=fontsize)

# Define tick positions to be at the center of each color segment
tick_positions = np.arange(len(terrain_types)) + 0.5  # Middle of each segment

# Adjust the colorbar boundaries to correspond to the color categories
boundaries = np.arange(len(terrain_types)+1)

# Add color bar with terrain names
cbar = plt.colorbar(sc, ax=ax, orientation='horizontal')
cbar.set_ticks(tick_positions)  
cbar.set_ticklabels(terrain_types, fontsize=fontsize-1.5)
cbar.set_label('Compositional Fields', fontsize=fontsize, labelpad=10)


# Make axes equal
ax.set_aspect('equal', adjustable='box')
ax.set_xlabel('x (km)', fontsize=fontsize)
ax.set_ylabel('y (km)', fontsize=fontsize)


# Save the plot
plt.savefig(plot_name, dpi=500)
plt.close()
