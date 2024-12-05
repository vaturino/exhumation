#! /usr/bin/python3
import numpy as np
import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from numpy.core.fromnumeric import size
from scipy.interpolate import griddata
from matplotlib.gridspec import GridSpec
import sys, os, subprocess
import json as json
from tqdm import tqdm
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'
import argparse
from pathlib import Path
import seaborn as sns
import matplotlib.tri as tri
path = str(Path(Path(__file__).parent.absolute()).parent.absolute())
sys.path.insert(0, path)
from libraries.functions import *
from libraries.particles import *
import plotly.graph_objects as go




def main():
    parser = argparse.ArgumentParser(description= 'Script that gets n models and the time at the end of computation and gives the temperature and viscosity plots for all time steps')
    parser.add_argument('json_file', help='json file with models and end time to plot T and eta field')  
    args = parser.parse_args()

    # csvs_loc =  '/home/vturino/PhD/projects/exhumation/gz_outputs/'
    # models_loc =  '/home/vturino/PhD/projects/exhumation/raw_outputs/'
    json_loc = '/home/vturino/PhD/projects/exhumation/pyInput/'

    models_loc = '/home/vturino/Vale_nas/exhumation/raw_outputs/'
    csvs_loc = '/home/vturino/Vale_nas/exhumation/gz_outputs/'

    with open(f"{json_loc}{args.json_file}") as json_file:
            configs = json.load(json_file)

    compositions = configs['compositions']
    cutoff = configs['cutoff']

    file_count = 0
    # Define the fixed color palette for the terrains
    terrain_palette = {
        "mantle": sns.color_palette("Paired")[0],
        "oc": "olivedrab",
        "sed": sns.color_palette("Paired")[2],
        "opc": sns.color_palette("Paired")[3],
        "serp": "darkslategray",
        "ecl": "sienna"  # Assign "sienna" color for "ecl"
    }

    # Define the fixed order of terrain types
    terrain_order = ["mantle", "oc", "sed", "opc", "serp", "ecl"]


    for ind_m, m in tqdm(enumerate(configs['models'])):    
        time_array = np.zeros((len(os.listdir(f"{csvs_loc}{m}/fields")),2)) 
        stat = pd.read_csv(f"{models_loc}{m}/statistics",skiprows=configs['head_lines'],sep='\s+',header=None)
        time_array = grab_dimTime_fields(f"{csvs_loc}{m}/fields", stat, time_array, configs['head_lines']-1)
        plot_loc_mod = f"/home/vturino/PhD/projects/exhumation/plots/single_models/{m}"
        if not os.path.exists(plot_loc_mod):
            os.mkdir(plot_loc_mod)
        plot_loc = f"{plot_loc_mod}/"
        if not os.path.exists(plot_loc):
            os.mkdir(plot_loc)
        else:
            # Count the files in the fields_loc directory
            file_count = len(os.listdir(plot_loc))

        
        for t in [0]:

            # Create the plot
            fig = plt.figure()
            gs = GridSpec(2, 1)
            plotname = f"{plot_loc}model_domain_with_contour.pdf"
            data = pd.read_parquet(f"{csvs_loc}{m}/fields/full.{int(t)}.gzip")
            data["lithology"] = 0
            data["Tcel"] = data["T"] - 273.15

            # Create lithology classes
            for ind, c in enumerate(compositions):
                data[c][data["Points:1"] < cutoff[ind]] = 0
                data[c][data[c] >= 0.5] = 1
                data[c][data[c] < 0.5] = 0

            for ind_c, c in enumerate(compositions):
                weight = ind_c + 1
                data["lithology"] += weight * data[c]

            data["lithology"] = data["lithology"].astype(int)
            composition_mapping = {ind_c + 1: c for ind_c, c in enumerate(compositions)}
            data["terrain"] = data["lithology"].map(composition_mapping)
            data["terrain"].fillna("mantle", inplace=True)

            xmin_plot = 0.
            xmax_plot = 5400.e3
            ymin_plot = 0.
            ymax_plot = 900.e3

            # Get present terrains in the current model
            present_terrains = [terrain for terrain in terrain_order if terrain in data["terrain"].unique()]

            # Create a color array in the fixed order
            colors_for_terrain = [terrain_palette[terrain] for terrain in terrain_order if terrain in present_terrains]

            # Map terrain values to numerical indices for coloring
            data["terrain_idx"] = data["terrain"].map({label: idx for idx, label in enumerate(present_terrains)})

            x = data["Points:0"].to_numpy() / 1.e3
            y = (ymax_plot - data["Points:1"].to_numpy()) / 1.e3

            triang = tri.Triangulation(x, y)
            colors = matplotlib.colormaps['Accent'].resampled(len(compositions) + 1)
            colors.set_bad('white')

            # Mask off unwanted triangles
            max_radius = 10
            triangles = triang.triangles
            xtri = x[triangles] - np.roll(x[triangles], 1, axis=1)
            ytri = y[triangles] - np.roll(y[triangles], 1, axis=1)
            maxi = np.max(np.sqrt(xtri**2 + ytri**2), axis=1)
            triang.set_mask(np.logical_and(maxi > max_radius, y[triangles][:, 1] < 90.))

            # Plot the terrain
            p1 = plt.tripcolor(triang, data["terrain_idx"], cmap=matplotlib.colors.ListedColormap(colors_for_terrain), shading='gouraud', vmax=len(compositions), vmin=0)

            # Add contour for logVisc at the specified level
            # plt.tricontour(triang, data["logVisc"], levels=[22], colors='red', linewidths=0.5, linestyles='solid')
            plt.tricontour(triang, data["Tcel"], levels=[1000], colors='navy', linewidths=0.2, linestyles='solid')

            # Add a colorbar for terrains
            cbar = plt.colorbar(p1, orientation='horizontal', aspect=50)
            cbar.set_ticks(np.arange(len(colors_for_terrain)))
            cbar.set_ticklabels(present_terrains)

            # Plot settings
            plt.ylim([(ymax_plot - ymin_plot) / 1.e3, -5])
            plt.xlim([xmin_plot / 1.e3, xmax_plot / 1.e3])
            plt.axhline(y=660, color='k', linestyle='--', linewidth=0.5)
            plt.gca().set_aspect('equal')

            #make boundary lines less thick
            thick = 0.3
            plt.gca().spines['top'].set_visible(False)
            plt.gca().spines['right'].set_linewidth(thick)
            plt.gca().spines['left'].set_linewidth(thick)
            plt.gca().spines['bottom'].set_linewidth(thick)
            plt.gca().tick_params(width=thick)

            # Save the plot
            plt.savefig(plotname, bbox_inches='tight', format='pdf', dpi=1000)
            plt.clf()
            plt.close('all')


if __name__ == "__main__":
    main()

