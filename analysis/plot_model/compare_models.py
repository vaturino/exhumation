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
from matplotlib.patches import Rectangle



def main():
    parser = argparse.ArgumentParser(description='Script that gets n models and the time at the end of computation and gives the temperature and viscosity plots for all time steps')
    parser.add_argument('json_file', help='json file with models and end time to plot T and eta field')  
    args = parser.parse_args()

    json_loc = '/home/vturino/PhD/projects/exhumation/pyInput/'
    models_loc = '/home/vturino/Vale_nas/exhumation/raw_outputs/'
    csvs_loc = '/home/vturino/Vale_nas/exhumation/gz_outputs/'

    with open(f"{json_loc}{args.json_file}") as json_file:
        configs = json.load(json_file)

    compositions = configs['compositions']
    cutoff = configs['cutoff']

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

    file_count = 0

    for ind_m, m in tqdm(enumerate(configs['models'])):
        time_array = np.zeros((len(os.listdir(f"{csvs_loc}{m}/fields")), 2)) 
        stat = pd.read_csv(f"{models_loc}{m}/statistics", skiprows=configs['head_lines'], sep='\s+', header=None)
        time_array = grab_dimTime_fields(f"{csvs_loc}{m}/fields", stat, time_array, configs['head_lines']-1)

        plot_loc_mod = f"/home/vturino/PhD/projects/exhumation/plots/single_models/{m}"
        if not os.path.exists(plot_loc_mod):
            os.mkdir(plot_loc_mod)
        plot_loc = f"{plot_loc_mod}/"
        if not os.path.exists(plot_loc):
            os.mkdir(plot_loc)
        else:
            file_count = len(os.listdir(plot_loc))

        f1, a1 = plt.subplots(3,3, figsize=(15, 9))
        plotname = f"{plot_loc}time_evolution.ps" 

    
        # Before the plotting loop, enable vector simplification globally
        plt.rcParams['path.simplify'] = True
        plt.rcParams['path.simplify_threshold'] = 0.5  # Tune as needed

        for indt, t in enumerate([1, 50, 90]):
            
            data = pd.read_parquet(f"{csvs_loc}{m}/fields/full.{int(t)}.gzip")
            data["lithology"] = 0
            data["logvisc"] = np.log10(data["viscosity"])
            data["logSR"] = np.log10(data["strain_rate"])

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

            filter_mask = (data["terrain"] != "mantle") & (data["terrain"] != "opc")
            data_filtered_comp = data[filter_mask]

            pts = get_points_with_y_in(data, 15.e3, 2.e3, ymax=900.e3)
            trench = get_trench_position(pts, threshold=0.13e7)
            xmin_plot = trench - 100.e3
            xmax_plot = trench + 200.e3
            ymin_plot = 740.e3
            ymax_plot = 902.e3

            x = data["Points:0"].to_numpy() / 1.e3
            y = (ymax_plot - data["Points:1"].to_numpy()) / 1.e3

            triang = tri.Triangulation(x, y)
            max_radius = 10
            triangles = triang.triangles

            xtri = x[triangles] - np.roll(x[triangles], 1, axis=1)
            ytri = y[triangles] - np.roll(y[triangles], 1, axis=1)
            maxi = np.max(np.sqrt(xtri**2 + ytri**2), axis=1)
            triang.set_mask(np.logical_and(maxi > max_radius, y[triangles][:, 1] < 90.))

            # Get present terrains in the current model
            present_terrains = [terrain for terrain in terrain_order if terrain in data["terrain"].unique()]

            # Create a color array in the fixed order
            colors_for_terrain = [terrain_palette[terrain] for terrain in terrain_order if terrain in present_terrains]

            # Map terrain values to numerical indices for coloring
            data["terrain_idx"] = data["terrain"].map({label: idx for idx, label in enumerate(present_terrains)})

            # Plot lithology using terrain indices and corresponding color map
            p1 = a1[0, indt].tripcolor(triang, data["terrain_idx"], cmap=matplotlib.colors.ListedColormap(colors_for_terrain), shading='gouraud', vmin=0, vmax=len(colors_for_terrain)-1, rasterized=True)
            a1[0, indt].tricontour(triang, data["T"] - 273.5, colors='k', levels=[100, 300, 500, 700, 900], linewidths=0.3)

            # Terrain colorbar


            a1[0, indt].set_ylim([(ymax_plot - ymin_plot) / 1.e3, -5])
            a1[0, indt].set_xlim([xmin_plot / 1.e3, xmax_plot / 1.e3])
            a1[0, indt].set_aspect('equal')
            a1[0, indt].spines[['top']].set_visible(False)
            # time_label = f"Time: {time_array[t, 1]/1e6:.1f} Myr"
            # a1[0, indt].text(0.05, 0.1, time_label, transform=a1[0, indt].transAxes, fontsize=12, verticalalignment='top', bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))



            # Strain rate plot
            p2 = a1[1, indt].tripcolor(triang, data["logSR"], cmap='RdBu_r', shading='gouraud', vmin=-19, vmax=-12, rasterized=True)
            a1[1, indt].spines[['top']].set_visible(False)
            a1[1, indt].set_ylim([(ymax_plot - ymin_plot) / 1.e3, -5])
            a1[1, indt].set_xlim([xmin_plot / 1.e3, xmax_plot / 1.e3])
            a1[1, indt].set_aspect('equal')
            # Annotate time in subplot a1[1, indt]
            

            # Viscosity plot + velocity vectors
            p3 = a1[2, indt].tripcolor(triang, data["logvisc"], shading='gouraud', vmin=18, vmax=24, rasterized=True)
            a1[2, indt].spines[['top']].set_visible(False)
            a1[2, indt].set_ylim([(ymax_plot - ymin_plot) / 1.e3, -5])
            a1[2, indt].set_xlim([xmin_plot / 1.e3, xmax_plot / 1.e3])
            a1[2, indt].set_aspect('equal')

        # Make all plots the same size
        f1.tight_layout()
        # Add axes labels
        a1[2,0].set_xlabel('x [km]')
        a1[2,1].set_xlabel('x [km]')
        a1[2,2].set_xlabel('x [km]')
        a1[0,0].set_ylabel('y [km]')
        a1[1,0].set_ylabel('y [km]')
        a1[2,0].set_ylabel('y [km]')



        # Colorbars
        cbar1 = f1.add_axes([0.92, 0.70, 0.02, 0.25])  # Adjust the position and size as needed
        cbar = plt.colorbar(p1, cax = cbar1, orientation='vertical')
        cbar.set_ticks(np.arange(len(colors_for_terrain)))
        cbar.set_ticklabels([terrain for terrain in present_terrains])  # Set tick labels to present terrain names
        cbar.set_label('Terrain', rotation=90, labelpad=0)  
        cbar2 = f1.add_axes([0.92, 0.385, 0.02, 0.25])
        plt.colorbar(p2, orientation='vertical', label='Log(Strain rate) [s-1]', cax=cbar2)
        cbar_ax = f1.add_axes([0.92, 0.07, 0.02, 0.25])  # Adjust the position and size as needed
        plt.colorbar(p3, cax=cbar_ax, orientation='vertical', label='Log(Viscosity) [Pa s]')
        plt.subplots_adjust(hspace = 0.1, right=0.9)  

        # # Enable vector simplification
        # plt.rcParams['path.simplify'] = True
        # plt.rcParams['path.simplify_threshold'] = 0.5  # Adjust threshold as needed


        plt.savefig(plotname, bbox_inches='tight', format='ps', dpi=500, pad_inches=0.1)
        plt.clf()
        plt.close('all')

if __name__ == "__main__":
    main()
