#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import pandas as pd
import argparse
import os, sys
import json
from pathlib import Path
from tqdm import tqdm
from matplotlib.colors import LinearSegmentedColormap
path = str(Path(Path(__file__).parent.absolute()).parent.absolute())
sys.path.insert(0, path)
from libraries.functions import *
from libraries.particles import *

def main():
    parser = argparse.ArgumentParser(description='Script to plot temperature and viscosity fields.')
    parser.add_argument('json_file', help='JSON file with models and end time to plot T and eta field')  
    args = parser.parse_args()

    json_loc = '/home/vturino/PhD/projects/exhumation/pyInput/'
    models_loc = '/home/vturino/Vale_nas/exhumation/raw_outputs/'
    csvs_loc = '/home/vturino/Vale_nas/exhumation/gz_outputs/'

    with open(f"{json_loc}{args.json_file}") as json_file:
        configs = json.load(json_file)

    # Define a single color (e.g., RGB format)
    single_color = ["gray", "forestgreen", "blue"]
    background_color = ["silver", "palegreen", "skyblue"]
    background_alpha = 0.2

   

    for ind_m, m in tqdm(enumerate(configs['models'])):    
        # Read in time step data
        file_count = len(os.listdir(f"{csvs_loc}{m}/fields"))
        plot_loc_mod = f"/home/vturino/PhD/projects/exhumation/plots/single_models/{m}"
        if not os.path.exists(plot_loc_mod):
            os.mkdir(plot_loc_mod)
        plot_loc = f"{plot_loc_mod}/Slab_profile/"
        if not os.path.exists(plot_loc):
            os.mkdir(plot_loc)

        fig, ax = plt.subplots(1, 3, figsize=(15, 5))
        plotname = f"{plot_loc}slab.pdf" 

        for indt, t in enumerate([1, 50, 90]):  # Loop over selected time steps
            tr = 1e-2
             # Create a colormap with just that color
            cmap = LinearSegmentedColormap.from_list("single_color", [single_color[indt], single_color[indt]])

            
            data = pd.read_parquet(f"{csvs_loc}{m}/fields/full.{int(t)}.gzip") 
            # data = data.sample(frac=0.1)  # Reduces the number of data points by half

            data["Tcel"] = data["T"] - 273.15  # Convert to Celsius
            #eliminate data["opc"] != 0 and data["op"] != 0
            # data = data[(data["opc"] < tr) & (data["op"] <tr)]


            pts = get_points_with_y_in(data, 15.e3, 2.e3, ymax=900.e3)
            trench = get_trench_position(pts, threshold=0.13e7)
            xmin_plot = trench - 400.e3
            xmax_plot = trench + 700.e3
            ymin_plot = 0.e3
            ymax_plot = 902.e3

            # Extract coordinates
            x = data["Points:0"].to_numpy() / 1.e3
            y = (ymax_plot - data["Points:1"].to_numpy()) / 1.e3
            triang = tri.Triangulation(x, y)

            # Mask off unwanted triangles based on max_radius
            max_radius = 10
            triangles = triang.triangles
            xtri = x[triangles] - np.roll(x[triangles], 1, axis=1)
            ytri = y[triangles] - np.roll(y[triangles], 1, axis=1)
            maxi = np.max(np.sqrt(xtri**2 + ytri**2), axis=1)
            triang.set_mask(np.logical_and(maxi > max_radius, y[triangles][:,1] < 90.))

            # if data Tcel > 700. put it to nan
            data["Tcel"] = np.where(data["Tcel"] > 1100, np.nan, data["Tcel"])

            # Mask values above 700°C (keep values < 700°C)
            Tcel_array = np.array(data["Tcel"])
           

            ax[indt].tripcolor(triang, Tcel_array, shading='flat', cmap=cmap, rasterized=True)
            ax[indt].axvspan(xmin_plot/1.e3, (xmax_plot-2000)/1.e3, color=background_color[indt], alpha=background_alpha)
            ax[indt].set_ylim(ymax_plot/1.e3, ymin_plot/1.e3)
            ax[indt].set_aspect('equal')
            ax[indt].axhline(y = (ymax_plot - 242.e3)/1e3, color = 'k', linestyle = '--', linewidth = 2)
            ax[indt].spines[['top']].set_visible(False)
            ax[indt].set_xticks([])
            ax[indt].set_yticks([])
            ax[indt].set_xlim(xmin_plot/1.e3, (xmax_plot-2000)/1.e3)

 


        fig.subplots_adjust(hspace=0.2)  
        fig.subplots_adjust(wspace=0.2)
        
        plt.savefig(plotname, bbox_inches='tight', format='pdf', dpi=500, transparent = True)
        plt.clf()
        plt.close('all')

if __name__ == "__main__":
    main()
