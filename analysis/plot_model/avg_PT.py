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
    part_marker = ["o", "s", "*", "D", "v", "^", "x", "p", "h", "H"]
    label_part = ["t1", "t2", "t3", "t4", "t5", "t6", "t7", "t8", "t9", "t10"]
    color_part = ["blue", "purple", "black"]

    file_count = 0

    for ind_m, m in tqdm(enumerate(configs['models'])):
        time_array = np.zeros((len(os.listdir(f"{csvs_loc}{m}/fields")), 2)) 
        stat = pd.read_csv(f"{models_loc}{m}/statistics", skiprows=configs['head_lines'], sep='\s+', header=None)
        time_array = grab_dimTime_fields(f"{csvs_loc}{m}/fields", stat, time_array, configs['head_lines']-1)

        plot_loc_mod = f"/home/vturino/PhD/projects/exhumation/plots/single_models/{m}/"
        if not os.path.exists(plot_loc_mod):
            os.mkdir(plot_loc_mod)
        plot_loc = f"{plot_loc_mod}/Compo_synthetic_parts/"
        if not os.path.exists(plot_loc):
            os.mkdir(plot_loc)
        else:
            file_count = len(os.listdir(plot_loc))

        cr = pd.read_csv(f"{plot_loc_mod}/txt_files/2D_v.txt", sep="\s+")
        cr["conv_rate"].iloc[0] = cr["conv_rate"].iloc[1]

        txt_loc = f"{plot_loc_mod}/txt_files/"
        pt_loc = f"{txt_loc}/PT"
        exhumed_list = pd.read_csv(f"{txt_loc}exhumed_particles_timing.txt", sep="\s+")
        exhumed_list_oc = exhumed_list[exhumed_list["lithology"] == "oc"].reset_index(drop=True)
  

        # Identify clumps of data based on 'tin' values with gaps in between
        clump_ranges = []
        sorted_tin_values = sorted(exhumed_list_oc['tin'].unique())
        
        # Iterate through sorted tin values to find gaps and define clumps
        start = sorted_tin_values[0]
        for i in range(1, len(sorted_tin_values)):
            if sorted_tin_values[i] - sorted_tin_values[i - 1] > 1:  # Define gap threshold as 1
                clump_ranges.append((start, sorted_tin_values[i - 1]))
                start = sorted_tin_values[i]
        clump_ranges.append((start, sorted_tin_values[-1]))  # Add the last clump
        


    exhumed_list_oc["time_interval"] = None
    exhumed_list_oc["clump"] = None

    # Assign time_interval and clump ID to particles
    for start, end in clump_ranges:
        mask = (exhumed_list_oc['tin'] >= start) & (exhumed_list_oc['tin'] <= end)
        exhumed_list_oc.loc[mask, 'time_interval'] = f"{start}-{end}"
        exhumed_list_oc.loc[mask, 'clump'] = clump_ranges.index((start, end)) + 1

    # Identify all unique clump IDs
    unique_clumps = exhumed_list_oc["clump"].dropna().unique()

    # Optional: Store synthetic particles for each clump
    all_synparts = {}

    # Loop over each clump
    for clump_id in unique_clumps:
        clump_particles = exhumed_list_oc[exhumed_list_oc["clump"] == clump_id]["id"].values
        combined_data = []

        # Collect data for all particles in this clump
        for p in clump_particles:
            pt_single = pd.read_csv(f"{pt_loc}/pt_part_{p}.txt", sep="\s+")
            combined_data.append(pt_single)

        # Concatenate all particle data and compute averages per time step
        combined_df = pd.concat(combined_data, axis=0, keys=clump_particles)
        avg_data = combined_df.groupby(level=1).mean()

        # Save synthetic particle data for this clump
        all_synparts[clump_id] = avg_data

        
    for clump_id, synpart in all_synparts.items():
        Pthresh = synpart["Plith"].max() * 0.75
        plt.plot(synpart["time"], synpart["Plith"], label=f"Burst {clump_id}", color=color_part[(clump_id - 1) % len(color_part)]) #, marker=part_marker[(clump_id - 1) % len(part_marker)], markersize=3)
        plt.axhline(y=Pthresh, color=color_part[(clump_id - 1) % len(color_part)], linestyle='--', linewidth=1, label =f"Threshold {clump_id}")
    plt.xlabel("Time (Myr)")
    plt.ylabel("Pressure (GPa)")
    # plt.grid()
    plt.legend()
    plt.title("Average Pressure of Synthetic Particles in Each Burst")
    plt.savefig(f"{plot_loc_mod}/average_pressure_synthetic_particles.pdf", dpi=500)
    plt.close()

    for t in tqdm(range(1, len(time_array))):
            f1, a1 = plt.subplots(2, 1, figsize=(8, 5), dpi=500, height_ratios=[1, 0.25])

            plotname = f"{plot_loc}{int(t)}.png" 
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

            present_terrains = [terrain for terrain in terrain_order if terrain in data["terrain"].unique()]
            colors_for_terrain = [terrain_palette[terrain] for terrain in terrain_order if terrain in present_terrains]
            data["terrain_idx"] = data["terrain"].map({label: idx for idx, label in enumerate(present_terrains)})

            # Plot lithology using terrain indices and corresponding color map
            p1 = a1[0].tripcolor(triang, data["terrain_idx"], cmap=matplotlib.colors.ListedColormap(colors_for_terrain), shading='gouraud', vmin=0, vmax=len(colors_for_terrain)-1)
            time_in_myr = time_array[t, 1] / 1.e6
            a1[0].annotate(f"Time = {time_in_myr:.1f} Myr", xy=(0.02, 0.05), xycoords='axes fraction', ha='left', fontsize=14, color='black')

            for clump_id, synpart in all_synparts.items():
                row = synpart.iloc[t-1]
                x_plot = row["x"] / 1.e3
                y_plot = (ymax_plot - row["y"]) / 1.e3
                a1[0].scatter(
                    x_plot,
                    y_plot,
                    s=20,
                    c=color_part[(clump_id - 1) % len(color_part)],
                    marker=part_marker[(clump_id - 1) % len(part_marker)],
                    zorder=100,
                    label=f"Burst {clump_id}"
                )

            matplotlib.rcParams['font.family'] = 'Arial'
            font_sz = 12
            cbar = plt.colorbar(p1, orientation='vertical', ax=a1[0], pad=0.02, shrink=0.95)
            cbar.set_ticks(np.arange(len(colors_for_terrain)))
            cbar.set_ticklabels([terrain for terrain in present_terrains])  # Set tick labels to present terrain names
            cbar.set_label('Lithology', rotation=90, labelpad=1, fontsize=font_sz)
            a1[0].set_ylim([(ymax_plot - ymin_plot) / 1.e3, -5])
            a1[0].set_xlim([xmin_plot / 1.e3, xmax_plot / 1.e3])
            a1[0].set_aspect('equal')
            a1[0].set_xlabel("x (km)", fontsize=font_sz)
            a1[0].set_ylabel("Depth (km)", fontsize=font_sz)
            a1[0].legend(loc='lower right', fontsize=font_sz, markerscale=1.5, frameon=False)
            a1[0].spines['top'].set_visible(False)
            a1[0].tick_params(axis='both', which='major', labelsize=font_sz)

            # plot convergence rate
            a1[1].plot(cr["time"]/1.e6, cr["conv_rate"], color="slategrey", linewidth=2, label="Convergence rate")
            a1[1].scatter(cr["time"].iloc[t]/1.e6, cr["conv_rate"].iloc[t], color="darkred", s=50, zorder=100, clip_on=False)
            # Add boxes around the 3 clumps in time
            y_min = 0 
            y_max = 8
            for ind, (clump_start, clump_end) in enumerate(clump_ranges):
                a1[1].add_patch(plt.Rectangle(
                    (clump_start, y_min),  # Bottom-left corner
                    (clump_end - clump_start),  # Width
                    y_max - y_min,  # Height
                    color=color_part[ind % len(color_part)], alpha=0.3, zorder=1  # Ensure valid index
                ))
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