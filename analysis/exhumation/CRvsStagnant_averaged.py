#!/usr/bin/env python3

import pandas as pd
import numpy as np
import json
import matplotlib.pyplot as plt
from pathlib import Path
import argparse
import os
from matplotlib.patches import Polygon
from sklearn.cluster import DBSCAN
from shapely.geometry import Polygon as ShapelyPolygon

############### MAIN ####################

def main():
    parser = argparse.ArgumentParser(description='Script to plot kinematic indicators with lithology hues.')
    parser.add_argument('json_file', help='JSON file with model name, time at the end of computation, folder where to save the plots, legend')
    args = parser.parse_args()

    json_loc = '/home/vturino/PhD/projects/exhumation/pyInput/'

    with open(f"{json_loc}{args.json_file}") as json_file:
        configs = json.load(json_file)
    m = configs["models"][0]

    plot_loc = f"/home/vturino/PhD/projects/exhumation/plots/single_models/{m}"
    txt_loc = f'{plot_loc}/txt_files'
    if not os.path.exists(txt_loc):
        os.mkdir(txt_loc)

    # Read the data from the CSV file
    data = pd.read_csv(f"{txt_loc}/stagnant_particles.txt", sep="\s+")
    cr = pd.read_csv(f"{txt_loc}/2D_v.txt", sep="\s+")
    cr.iloc[0] = np.nan
    cr = cr.dropna()

    # Define colors for lithologies
    colors_tmax = {
        "sed": "mediumblue",
        "oc": "#B06D1A",
        "ecl": "#45701C",
        "serp": "brown"
    }

    colors_tfin = {
        "sed": "cornflowerblue",
        "oc": "#E3B64F",
        "ecl": "#A0C93D",
        "serp": "lightsalmon"
    }

    # Initialize storage for extremes
    extremes_dict = {"kin": {}, "dyn": {}, "trans": {}}

    # Plot
    f1, a1 = plt.subplots(2, 1, figsize=(15, 7), height_ratios=[0.25, 1])

    # Plot convergence rate in the first subplot
    a1[0].plot(cr['time'] / 1e6, cr['conv_rate'], color='grey', label='Convergence rate', linewidth=2)
    a1[0].label_outer()
    a1[0].set_xlim(0, 50)
    a1[0].set_ylim(0, max(cr['conv_rate']) + 0.2)
    a1[0].axvline(x=35., color='grey', linewidth=1, linestyle="--")
    a1[0].set_title('Convergence Rate')

    # Plot stagnation intervals and scatter points in the second subplot
    for lithology in data['lithology'].unique():
        subset = data[data['lithology'] == lithology]

        # Scatter points for tm_kin, tm_dyn, and tm_trans
        a1[1].scatter(
            subset['tm_kin'], subset['Pm_kin'],
            color=colors_tmax[lithology], label=f'{lithology} (kin)', alpha=1, s=5, zorder=10
        )
        a1[1].scatter(
            subset['tm_dyn'], subset['Pm_dyn'],
            color=colors_tmax[lithology], label=f'{lithology} (dyn)', alpha=1, s=5, zorder=10
        )
        a1[1].scatter(
            subset['tm_trans'], subset['Pm_trans'],
            color=colors_tmax[lithology], label=f'{lithology} (trans)', alpha=1, s=5, zorder=10
        )

        # Collect extreme points for each time stage and store them in dictionary
        for _, row in subset.iterrows():
            # Collect extreme points for each time stage
            kin_extremes = [
                (row['ti_kin'], row['Pm_kin']),
                (row['ti_kin'] + row['time_interval_kin'], row['Pm_kin']),
            ]
            dyn_extremes = [
                (row['ti_dyn'], row['Pm_dyn']),
                (row['ti_dyn'] + row['time_interval_dyn'], row['Pm_dyn']),
            ]
            trans_extremes = [
                (row['ti_trans'], row['Pm_trans']),
                (row['ti_trans'] + row['time_interval_trans'], row['Pm_trans']),
            ]
            
            # Store the extremes for later use in a dictionary by lithology and time stage
            if lithology not in extremes_dict["kin"]:
                extremes_dict["kin"][lithology] = []
            if lithology not in extremes_dict["dyn"]:
                extremes_dict["dyn"][lithology] = []
            if lithology not in extremes_dict["trans"]:
                extremes_dict["trans"][lithology] = []

            extremes_dict["kin"][lithology].extend(kin_extremes)
            extremes_dict["dyn"][lithology].extend(dyn_extremes)
            extremes_dict["trans"][lithology].extend(trans_extremes)

    # Now plot the polygons for each lithology and each time stage (kin, dyn, trans)
    for time_stage in ['kin', 'dyn', 'trans']:
        for lithology, extremes in extremes_dict[time_stage].items():
            # Sort the extremes by pressure (ascending order)
            extremes_sorted = sorted(extremes, key=lambda x: x[1])  # Sort by pressure

            # Remove NaN values from pressures and times (if any)
            cleaned_extremes = [(time, pressure) for time, pressure in extremes_sorted if not np.isnan(time) and not np.isnan(pressure)]

            # If no valid points, skip this cluster
            if len(cleaned_extremes) < 4:
                continue

            # Separate times and pressures for the polygon
            polygon_x = [ext[0] for ext in cleaned_extremes]
            polygon_y = [ext[1] for ext in cleaned_extremes]

            # Clustering pressures using DBSCAN
            db = DBSCAN(eps=0.1, min_samples=3).fit(np.array(polygon_y).reshape(-1, 1))
            clusters = db.labels_

            # Plot polygons for each cluster
            for cluster_label in np.unique(clusters):
                if cluster_label == -1:
                    # -1 means noise, so we ignore it
                    continue

                # Collect points belonging to this cluster
                cluster_points = [cleaned_extremes[i] for i in range(len(clusters)) if clusters[i] == cluster_label]

                # Only proceed if the cluster has at least 4 points
                if len(cluster_points) < 4:
                    continue

                # Sort the points by pressure (ascending order)
                cluster_points_sorted = sorted(cluster_points, key=lambda x: x[1])

                # Separate times and pressures for the polygon
                polygon_x = [point[0] for point in cluster_points_sorted]
                polygon_y = [point[1] for point in cluster_points_sorted]

                # Create a Shapely polygon
                shapely_polygon = ShapelyPolygon(zip(polygon_x, polygon_y))

                # Simplify the polygon to reduce vertices while maintaining the shape
                simplified_polygon = shapely_polygon.simplify(0.05, preserve_topology=True)

                # Add the simplified polygon to the plot
                a1[1].add_patch(Polygon(list(simplified_polygon.exterior.coords), closed=True,
                                        color=colors_tfin[lithology], alpha=0.5, label=f'{lithology} ({time_stage})'))

    # Customize plot
    a1[1].set_xlabel('Time (Myr)')
    a1[1].set_ylabel('Pressure (GPa)')
    a1[1].set_ylim(0, 2.6)
    a1[1].set_xlim(0, 50)
    a1[1].axvline(x=35., color='grey', linewidth=1, linestyle="--")
    a1[1].set_title('Stagnation Time Range by Lithology')
    a1[1].legend(loc='upper right', fontsize=8, frameon=False)

    # Adjust spacing
    plt.subplots_adjust(hspace=0.1)

    # Save plot
    plt.savefig(f"{plot_loc}/stagnant_time_intervals_by_lithology_separate_times.png")
    plt.close()

if __name__ == "__main__":
    main()
