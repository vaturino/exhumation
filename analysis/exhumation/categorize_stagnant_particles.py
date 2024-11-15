#! /usr/bin/python3

import pandas as pd
import numpy as np
import json
import matplotlib.pyplot as plt
from pathlib import Path
import argparse
import os

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

        # Plot horizontal lines representing the time intervals
        for _, row in subset.iterrows():
            a1[1].plot(
                [row['ti_kin'], row['ti_kin'] + row['time_interval_kin']],
                [row['Pm_kin'], row['Pm_kin']],
                color=colors_tfin[lithology], linewidth=0.8, alpha=0.3
            )
            a1[1].plot(
                [row['ti_dyn'], row['ti_dyn'] + row['time_interval_dyn']],
                [row['Pm_dyn'], row['Pm_dyn']],
                color=colors_tfin[lithology], linewidth=0.8, alpha=0.3
            )
            a1[1].plot(
                [row['ti_trans'], row['ti_trans'] + row['time_interval_trans']],
                [row['Pm_trans'], row['Pm_trans']],
                color=colors_tfin[lithology], linewidth=0.8, alpha=0.3
            )

        # Scatter points for tm_kin, tm_dyn, and tm_trans
        a1[1].scatter(
            subset['tm_kin'], subset['Pm_kin'],
            color=colors_tmax[lithology], label=f'{lithology} (kin)', alpha=1, s=5, zorder = 10
        )
        a1[1].scatter(
            subset['tm_dyn'], subset['Pm_dyn'],
            color=colors_tmax[lithology], label=f'{lithology} (dyn)', alpha=1, s=5, zorder = 10
        )
        a1[1].scatter(
            subset['tm_trans'], subset['Pm_trans'],
            color=colors_tmax[lithology], label=f'{lithology} (trans)', alpha=1, s=5,  zorder = 10
        )

    # Customize plot
    a1[1].set_xlabel('Time (Myr)')
    a1[1].set_ylabel('Pressure (GPa)')
    a1[1].set_ylim(0, 2.6)
    # a1[1].axhline(y=0., color='grey', linewidth=1, linestyle="--")
    a1[1].set_xlim(0, 50)
    a1[1].axvline(x=35., color='grey', linewidth=1, linestyle="--")
    a1[1].set_title('Stagnation Time Range by Lithology')
    a1[1].legend(loc='upper right', fontsize=8, frameon=False)

    # Adjust spacing
    plt.subplots_adjust(hspace=0.1)

    # Save plot
    plt.savefig(f"{plot_loc}/stagnant_time_intervals_by_lithology.png")
    plt.close()


    # The lines in this plot are too many. I would like to plot, for each lithology, 3 lines: 
    # I have 3 values of each pressure: "shallow", intermediate" and "deep".
    # I would like to plot the avg time interval for each of these 3 values of pressure, for each lithology and each dyn, kin, trans.



if __name__ == "__main__":
    main()
