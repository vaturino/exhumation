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
#suppress pandas warnings
import warnings
warnings.filterwarnings("ignore")


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

    time_array = np.zeros((101, 2))
    time_array[:, 0] = np.linspace(0, 100, 101)  
    time_array[:, 1] = np.linspace(0, 50, 101) 

    plot_loc = f"/home/vturino/PhD/projects/exhumation/plots/comparison/{configs['test']}/"
    if not os.path.exists(plot_loc):
        os.makedirs(plot_loc)


    for t in tqdm(range(0, 30)):
        f1, a1 = plt.subplots(2, 3, figsize=(20, 10))
        plotname = f"{plot_loc}{int(t)}.png" 
        
        for ind_m, m in enumerate(configs['models']):
            data = pd.read_parquet(f"{csvs_loc}{m}/fields/full.{int(t)}.gzip")
            data["lithology"] = 0
            data["logSR"] = np.log10(data["strain_rate"])
            data["comp"] = data["oc"]+data["sed"]+data["ecl"]

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

            present_terrains = [terrain for terrain in terrain_order if terrain in data["terrain"].unique()]
            colors_for_terrain = [terrain_palette[terrain] for terrain in terrain_order if terrain in present_terrains]
            data["terrain_idx"] = data["terrain"].map({label: idx for idx, label in enumerate(present_terrains)})

            # Plot lithology using terrain indices and corresponding color map
            p1 = a1[0,ind_m].tripcolor(triang, data["terrain_idx"], cmap=matplotlib.colors.ListedColormap(colors_for_terrain), shading='gouraud', vmin=0, vmax=len(colors_for_terrain)-1)
            a1[0,ind_m].spines[['top']].set_visible(False)
        
            cbar = plt.colorbar(p1, orientation='horizontal', ax=a1[0,ind_m])
            cbar.set_ticks(np.arange(len(colors_for_terrain)))
            cbar.set_ticklabels([terrain for terrain in present_terrains])  # Set tick labels to present terrain names
            cbar.set_label('Terrain', rotation=0, labelpad=5)

            a1[0,ind_m].set_ylim([(ymax_plot - ymin_plot) / 1.e3, -5])
            a1[0,ind_m].set_xlim([xmin_plot / 1.e3, xmax_plot / 1.e3])
            a1[0,ind_m].set_aspect('equal')
            a1[0,ind_m].text(0.05, 0.1, configs['legend'][ind_m], transform=a1[0,ind_m].transAxes, fontsize=12, verticalalignment='top', bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))

            # Strain rate plot
            p2 = a1[1,ind_m].tripcolor(triang, data["logSR"], cmap='RdBu_r', shading='gouraud', vmin=-19, vmax=-12)
            a1[1,ind_m].tricontour(triang, data["comp"], colors='k', levels=[1], linewidths=0.7, alpha = 0.5)
            a1[1,ind_m].spines[['top']].set_visible(False)
            plt.colorbar(p2, orientation='horizontal', label='Log(Strain rate) [s-1]', ax=a1[1,ind_m])
            a1[1,ind_m].set_ylim([(ymax_plot - ymin_plot) / 1.e3, -5])
            a1[1,ind_m].set_xlim([xmin_plot / 1.e3, xmax_plot / 1.e3])
            a1[1,ind_m].set_aspect('equal')
            # Annotate time in subplot a1[1,ind_m]
            time_label = f"Time: {time_array[t, 1]:.1f} Myr"
            a1[1,ind_m].text(0.05, 0.1, time_label, transform=a1[1,ind_m].transAxes, fontsize=12, verticalalignment='top', bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))

        
        # Save and close the plot
        # print(f"Saving plot to {plotname}")
        plt.savefig(plotname, bbox_inches='tight', format='png', dpi=500)
        plt.clf()
        plt.close('all')

if __name__ == "__main__":
    main()
