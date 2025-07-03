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

        plot_loc_mod = f"/home/vturino/PhD/projects/exhumation/plots/single_models/{m}/"
        if not os.path.exists(plot_loc_mod):
            os.mkdir(plot_loc_mod)
        plot_loc = f"{plot_loc_mod}/Compo/"
        if not os.path.exists(plot_loc):
            os.mkdir(plot_loc)
        else:
            file_count = len(os.listdir(plot_loc))

        cr = pd.read_csv(f"{plot_loc_mod}/txt_files/2D_v.txt", sep="\s+")
        cr["conv_rate"].iloc[0] = cr["conv_rate"].iloc[1]

        for t in tqdm(range(file_count, len(time_array))):
        # for t in tqdm(range(60, 80)):
            f1, a1 = plt.subplots(2, 1, figsize=(8, 5), dpi=500, height_ratios=[1, 0.25])

            plotname = f"{plot_loc}{int(t)}.png" 
            # plotname = f"{plot_loc}{t}.pdf"
            data = pd.read_parquet(f"{csvs_loc}{m}/fields/full.{int(t)}.gzip")
            data["lithology"] = 0
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
            data_filtered_comp = data[filter_mask]

            pts = get_points_with_y_in(data, 15.e3, 2.e3, ymax=900.e3)
            trench = get_trench_position(pts, threshold=0.13e7)
            xmin_plot = trench - 100.e3
            xmax_plot = trench + 200.e3
            ymin_plot = 760.e3
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
            
            # annotate time in Myr with 1 decimal
            time_in_myr = time_array[t, 1] / 1.e6
            a1[0].annotate(f"Time = {time_in_myr:.1f} Myr", xy=(0.02, 0.05), xycoords='axes fraction', ha='left', fontsize=14, color='black')

            # Set Calibri as the font
            matplotlib.rcParams['font.family'] = 'Arial'
            font_sz = 12



           
            a1[0].clabel(
                a1[0].tricontour(triang, data["T"] - 273.5, colors='k', levels=[300, 700, 1100], linewidths=0.5),
                inline=False,
                fontsize=10,
                fmt='%1.0f°C',
                manual = [(trench / 1.e3 + 170, 0),   (trench / 1.e3 + 170, 35),   (trench / 1.e3 + 170, 70)],

                use_clabeltext=True
            )
            for label in a1[0].texts:
                label.set_rotation(0)
                label.set_fontweight('bold')
                label.set_bbox(dict(facecolor='white', edgecolor='none', alpha=0.5, boxstyle='round,pad=0.2'))

     




            # Terrain colorbar (vertical on the right)
            cbar = plt.colorbar(p1, orientation='vertical', ax=a1[0], pad=0.02, shrink=0.95)
            cbar.set_ticks(np.arange(len(colors_for_terrain)))
            cbar.set_ticklabels([terrain for terrain in present_terrains])  # Set tick labels to present terrain names
            cbar.set_label('Lithology', rotation=90, labelpad=1, fontsize=font_sz)
            a1[0].set_ylim([(ymax_plot - ymin_plot) / 1.e3, -5])
            a1[0].set_xlim([xmin_plot / 1.e3, xmax_plot / 1.e3])
            a1[0].set_aspect('equal')
            a1[0].set_xlabel("x (km)", fontsize=font_sz)
            a1[0].set_ylabel("Depth (km)", fontsize=font_sz)

            # Remove the top spine
            a1[0].spines['top'].set_visible(False)

            # Set font size for tick labels
            a1[0].tick_params(axis='both', which='major', labelsize=font_sz)

            # plot convergence rate
            a1[1].plot(cr["time"]/1.e6, cr["conv_rate"], color="slategrey", linewidth=2, label="Convergence rate")
            a1[1].scatter(cr["time"].iloc[t]/1.e6, cr["conv_rate"].iloc[t], color="darkred", s=50, zorder=100, clip_on=False)
            a1[1].set_xlim([0, 52])
            a1[1].set_ylim([0, 8])
            a1[1].set_xlabel("Time (Myr)", fontsize=font_sz)
            a1[1].set_ylabel(r"$v_c$ (cm/yr)", fontsize=font_sz)

            # Set font size for tick labels
            a1[1].tick_params(axis='both', which='major', labelsize=font_sz)

            # Hide all spines
            for side in ["top", "right", "bottom", "left"]:
                a1[1].spines[side].set_visible(False)

            # Hide tick marks but keep tick labels
            a1[1].tick_params(left=True, bottom=True)

            # Get axis limits to position arrows properly
            x_min, x_max = a1[1].get_xlim() 
            y_min, y_max = a1[1].get_ylim()

            x_min = x_min - 0.5
            y_min = y_min - 0.8

            # Add custom x-axis with arrow
            a1[1].annotate('', xy=(x_max, 0), xytext=(x_min, 0),
                        arrowprops=dict(arrowstyle='->', linewidth=1),
                        xycoords='data', zorder=5)

            # Add custom y-axis with arrow
            a1[1].annotate('', xy=(0, y_max), xytext=(0, y_min),
                        arrowprops=dict(arrowstyle='->', linewidth=1),
                        xycoords='data', zorder=5)
            

            

            # Save and close the plot
            plt. subplots_adjust(hspace=0.1)
            plt.tight_layout()
            plt.savefig(plotname, format='png', dpi=500)
            plt.clf()
            plt.close('all')

if __name__ == "__main__":
    main()
