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

pd.options.mode.chained_assignment = None  # Suppress warnings about chained assignments


def main():
    # Parse input arguments
    parser = argparse.ArgumentParser(description="Generate strain rate plots with different colors for plastic_yielding values.")
    parser.add_argument('json_file', help='JSON file containing model configurations.')
    args = parser.parse_args()

    # Define directories
    json_loc = '/home/vturino/PhD/projects/exhumation/pyInput/'
    models_loc = '/home/vturino/Vale_nas/exhumation/raw_outputs/'
    csvs_loc = '/home/vturino/Vale_nas/exhumation/gz_outputs/'

    # Load JSON configuration
    with open(f"{json_loc}{args.json_file}") as json_file:
        configs = json.load(json_file)

    # Process each model
    for ind_m, model in tqdm(enumerate(configs['models'])):
        plot_loc = f"/home/vturino/PhD/projects/exhumation/plots/single_models/{model}/Plastic_vs_Viscous/"
        os.makedirs(plot_loc, exist_ok=True)

        # Determine timesteps
        time_array = np.zeros((len(os.listdir(f"{csvs_loc}{model}/fields")), 2))

        for t in tqdm(range(0, len(time_array), 2)):
            fig, ax = plt.subplots()
            plotname = f"{plot_loc}{int(t/2)}.png"
            data = pd.read_parquet(f"{csvs_loc}{model}/fields/full.{int(t)}.gzip")
            data["logSR"] = np.log10(data["strain_rate"])
            
            pts = get_points_with_y_in(data, 15.e3, 2.e3, ymax=900.e3)
            trench = get_trench_position(pts, threshold=0.13e7)
            xmin_plot = trench - 100.e3
            xmax_plot = trench + 200.e3
            ymin_plot = 740.e3
            ymax_plot = 900.e3

            x = data["Points:0"].to_numpy() / 1.e3
            y = (ymax_plot - data["Points:1"].to_numpy()) / 1.e3

            # Define max_radius for triangle side length filtering
            max_radius = 10

            # Separate masks for plastic_yielding values
            mask_0 = data["plastic_yielding"] == 0
            mask_1 = data["plastic_yielding"] == 1


           # Generate and mask triangulations for plastic_yielding == 0
            if mask_0.sum() > 0:
                triang_0 = tri.Triangulation(x[mask_0], y[mask_0])
                xtri_0 = x[mask_0][triang_0.triangles] - np.roll(x[mask_0][triang_0.triangles], 1, axis=1)
                ytri_0 = y[mask_0][triang_0.triangles] - np.roll(y[mask_0][triang_0.triangles], 1, axis=1)
                side_lengths_0 = np.sqrt(xtri_0**2 + ytri_0**2)
                max_side_0 = np.max(side_lengths_0, axis=1)
                triang_0.set_mask(max_side_0 > max_radius)

            # Generate and mask triangulations for plastic_yielding == 1
            if mask_1.sum() > 0:
                triang_1 = tri.Triangulation(x[mask_1], y[mask_1])
                xtri_1 = x[mask_1][triang_1.triangles] - np.roll(x[mask_1][triang_1.triangles], 1, axis=1)
                ytri_1 = y[mask_1][triang_1.triangles] - np.roll(y[mask_1][triang_1.triangles], 1, axis=1)
                side_lengths_1 = np.sqrt(xtri_1**2 + ytri_1**2)
                max_side_1 = np.max(side_lengths_1, axis=1)
                triang_1.set_mask(max_side_1 > max_radius)

            # Plot strain_rate for plastic_yielding == 0
            if mask_0.sum() > 0:
                ax.tripcolor(triang_0, data["logSR"][mask_0], cmap="Reds", shading="gouraud", vmin=-19, vmax=-12)

            # Plot strain_rate for plastic_yielding == 1
            if mask_1.sum() > 0:
                ax.tripcolor(triang_1, data["logSR"][mask_1], cmap="Blues", shading="gouraud", vmin=-19, vmax=-12)

                # Add contour for plastic_yielding == 1
                if mask_1.sum() > 0:  # Check if there are any points where plastic_yielding == 1
                    # Plot the contour for plastic_yielding == 1
                    contour = ax.tricontour(
                        triang_1, 
                        data["plastic_yielding"][mask_1],  # Use the same mask for contouring
                        levels=[1],  # Only show the contour where plastic_yielding == 1
                        colors='black',  # Black color for the contour
                        linewidths=1.5  # Adjust line width for visibility
                    )

            # Set plot limits and aspect ratio
            ax.set_ylim([(ymax_plot - ymin_plot) / 1.e3, -5])
            ax.set_xlim([xmin_plot / 1.e3, xmax_plot / 1.e3])
            ax.set_aspect('equal')
            ax.spines['top'].set_visible(False)

            # Save the plot
            plt.savefig(plotname, bbox_inches='tight', format='png', dpi=500)
            plt.clf()
            plt.close('all')


if __name__ == "__main__":
    main()
