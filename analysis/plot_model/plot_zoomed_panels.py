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
        plot_loc = f"{plot_loc_mod}/Panels_pdf/"
        if not os.path.exists(plot_loc):
            os.mkdir(plot_loc)
        else:
            file_count = len(os.listdir(plot_loc))

        # for t in tqdm(range(0, len(time_array), 2)):
        for t in [50, 70, 90]:
            f1, a1 = plt.subplots(1, 3, figsize=(20, 5))
            plotname = f"{plot_loc}{int(t/2)}.pdf" 
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
            p1 = a1[0].tripcolor(triang, data["terrain_idx"], cmap=matplotlib.colors.ListedColormap(colors_for_terrain), shading='gouraud', vmin=0, vmax=len(colors_for_terrain)-1)
            a1[0].tricontour(triang, data["T"] - 273.5, colors='k', levels=[100, 300, 500, 700, 900], linewidths=0.3)

            # Terrain colorbar
            cbar = plt.colorbar(p1, orientation='horizontal', ax=a1[0])
            cbar.set_ticks(np.arange(len(colors_for_terrain)))
            cbar.set_ticklabels([terrain for terrain in present_terrains])  # Set tick labels to present terrain names
            cbar.set_label('Terrain', rotation=0, labelpad=5)

            a1[0].set_ylim([(ymax_plot - ymin_plot) / 1.e3, -5])
            a1[0].set_xlim([xmin_plot / 1.e3, xmax_plot / 1.e3])
            a1[0].set_aspect('equal')

            # Strain rate plot
            p2 = a1[1].tripcolor(triang, data["logSR"], cmap='RdBu_r', shading='gouraud', vmin=-19, vmax=-12)
            a1[1].spines[['top']].set_visible(False)
            plt.colorbar(p2, orientation='horizontal', label='Log(Strain rate) [s-1]', ax=a1[1])
            a1[1].set_ylim([(ymax_plot - ymin_plot) / 1.e3, -5])
            a1[1].set_xlim([xmin_plot / 1.e3, xmax_plot / 1.e3])
            a1[1].set_aspect('equal')
            # Annotate time in subplot a1[1]
            time_label = f"Time: {time_array[t, 1]/1e6:.1f} Myr"
            a1[1].text(0.05, 0.1, time_label, transform=a1[1].transAxes, fontsize=12, verticalalignment='top', bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))

            # Viscosity plot + velocity vectors
            p3 = a1[2].tripcolor(triang, data["logvisc"], shading='gouraud', vmin=18, vmax=24)
            vel_plot_thresh = 0.0001  # don't plot velocity vectors smaller than this (cm/yr)
            step = 500  # plot every 100th vector to reduce the number of vectors plotted
            # # Filter out small velocity vectors and downsample data for plotting
            # mask = (100. * np.sqrt(data["velocity:0"]**2 + data["velocity:1"]**2)) >= vel_plot_thresh
            # mask = mask.reindex(data_filtered_comp.index, fill_value=False)
            # data_filtered = data_filtered_comp[mask].iloc[::step]

            # vel_vects = a1[2].quiver(data_filtered["Points:0"].to_numpy() / 1.e3, 
            #                          (ymax_plot - data_filtered["Points:1"].to_numpy()) / 1.e3, 
            #                          data_filtered["velocity:0"].to_numpy() * 100, 
            #                          data_filtered["velocity:1"].to_numpy() * 100, 
            #                          scale=50, color='black', width=0.0015)
            # a1[2].quiverkey(vel_vects, 0.15, 0.1, 1, '1 cm/yr', labelpos='W', fontproperties={'size': '7'}, color='white', labelcolor='white')
            a1[2].spines[['top']].set_visible(False)
            plt.colorbar(p3, orientation='horizontal', label='Log(Viscosity) [Pa s]', ax=a1[2])
            a1[2].set_ylim([(ymax_plot - ymin_plot) / 1.e3, -5])
            a1[2].set_xlim([xmin_plot / 1.e3, xmax_plot / 1.e3])
            a1[2].set_aspect('equal')

            # Save and close the plot
            plt.savefig(plotname, bbox_inches='tight', format='pdf', dpi=500)
            plt.clf()
            plt.close('all')

if __name__ == "__main__":
    main()
