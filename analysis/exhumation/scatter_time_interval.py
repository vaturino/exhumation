#!/usr/bin/env python3

import pandas as pd
import numpy as np
import json
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import os
from matplotlib.colors import Normalize
from matplotlib import cm
from matplotlib import gridspec

def main():
    # Parse the arguments
    parser = argparse.ArgumentParser(description='Script that gets some models and gives the kinematic indicators')
    parser.add_argument('json_file', help='json file with model name, time at the end of computation, folder where to save the plots, legend')
    args = parser.parse_args()

    colors_tfin = {
    "sed": "cornflowerblue",
    "oc": "#E3B64F",
    "ecl": "#A0C93D",
    "serp": "lightsalmon"
    }

    json_loc = '/home/vturino/PhD/projects/exhumation/pyInput/'
    json_file = f"{json_loc}{args.json_file}"

    # Read the json file
    with open(f"{json_file}") as json_file:
        configs = json.load(json_file)

    # Create the folders to save the plots
    plot_loc = f"/home/vturino/PhD/projects/exhumation/plots/single_models/{configs['models'][0]}"
    txt_loc = f'{plot_loc}/txt_files'
    if not os.path.exists(txt_loc):
        os.mkdir(txt_loc)

    # Load the data
    exhumed_list = pd.read_csv(f"{txt_loc}/timing_exhumed_particles.txt", sep="\s+")
    stagnant_list = pd.read_csv(f"{txt_loc}/timing_stagnant_particles.txt", sep="\s+")
    cr = pd.read_csv(f"{txt_loc}/2D_v.txt", sep="\s+")
    cr["conv_rate"].iloc[0] = np.nan

    # Load additional timing data
    exhumed_list = pd.read_csv(f"{txt_loc}/exhumed_times.txt", sep="\s+")
    stagnant_list = pd.read_csv(f"{txt_loc}/stagnant_times.txt", sep="\s+")

    # Select plot parameters

    linewidth = 1
    rugline = 3
    rugalpha = 0.3

    # Create the figure and gridspec for subplots with space for colorbars
    fig,ax = plt.subplots(3, 1, figsize=(15, 10), height_ratios=[0.25, 1, 1])


    # Plot with the time interval duration for exhumed and stagnant vs conv rate
    ax[0].plot(cr["time"]/1.e6, cr["conv_rate"], color="grey")
    ax[0].set_ylim(0, 7.5)
    ax[0].set_xlim(0, 50)
    ax[0].axvline(x=35, color="grey", linestyle="--", linewidth=linewidth)
    ax[0].set_xticks([])
    ax[0].set_yticks([1, 3, 5, 7])

    # Plot Exhumed data
    sns.scatterplot(ax=ax[1], data=exhumed_list, x="ti", y="time_interval", hue="lithology", palette=colors_tfin, zorder=10, 
                    style="lithology", s = 100, markers = ["o", "d", "s", "X"], linewidth=0.15, edgecolor="black", alpha=0.8)
    sns.rugplot(ax=ax[1], data=exhumed_list, x="ti", y="time_interval", hue = "lithology", palette=colors_tfin, linewidth=rugline, alpha=rugalpha)
    # ax1 = ax[1].twiny()
    # sns.histplot(ax=ax1, data=exhumed_list, y="time_interval", hue = "lithology", palette=colors_tfin, bins=50, alpha=0.3, kde=False, zorder = 1)
    sns.kdeplot(ax=ax[1], data=exhumed_list, x="ti", y="time_interval", levels=10, color="magenta", linewidths=linewidth, alpha=0.5, fill=True)
    # sns.regplot(ax=ax[1], data=exhumed_list, x="ti", y="time_interval", scatter=False, color="magenta", ci=None, line_kws={"linewidth": linewidth}, robust=True)
    ax[1].set_xlim(0, 50)
    ax[1].axvline(x=35, color="grey", linestyle="--", linewidth=linewidth)
    ax[1].set_xticks([])
    ax[1].set_ylabel("Duration (Ma)")
    ax[1].set_xlabel("")
    ax[1].set_ylim(0, 35)
    ax[1].set_yticks([5, 10, 15, 20, 25, 30, 35])

    # Plot Stagnant data (Dynamic slowdown)
    # sns.scatterplot(ax=ax[2], data=stagnant_list, x="ti_dyn", y="time_interval_dyn", hue="lithology", palette=colors_tfin, zorder=10, 
    #                 style="lithology", s = 80, markers = ["o", "d", "s", "X"], linewidth=0.15, edgecolor="black", alpha=0.5)
    sns.kdeplot(ax=ax[2], data=stagnant_list, x="ti_dyn", y="time_interval_dyn", levels=10, color="grey", linewidths=linewidth, alpha=0.5, fill=True)   
    sns.rugplot(ax=ax[2], data=stagnant_list, x="ti_dyn", y="time_interval_dyn", hue = "lithology", palette=colors_tfin, linewidth=rugline, alpha=rugalpha) 
    
    # Plot Stagnant data (Kinematic slowdown)
    sns.scatterplot(ax=ax[2], data=stagnant_list, x="ti_kin", y="time_interval_kin", hue="lithology", palette=colors_tfin, zorder=10, 
                    style="lithology", s = 80, markers = ["o", "d", "s", "X"], linewidth=0.15, edgecolor="black", alpha=0.5)
    sns.kdeplot(ax=ax[2], data=stagnant_list, x="ti_kin", y="time_interval_kin", levels=10, color="indigo", linewidths=linewidth, alpha=0.5, fill=True)
    sns.rugplot(ax=ax[2], data=stagnant_list, x="ti_kin", y="time_interval_kin", hue = "lithology", palette=colors_tfin, linewidth=rugline, alpha=rugalpha)
    # sns.regplot(ax=ax[2], data=stagnant_list, x="ti_kin", y="time_interval_kin", scatter=False, color="blue", line_kws={"linewidth": linewidth})

    # Plot transient point data
    sns.scatterplot(ax=ax[2], data=stagnant_list, x="ti_trans", y="time_interval_trans", hue="lithology", palette=colors_tfin, zorder=10,
                    style="lithology", s = 80, markers = ["o", "d", "s", "X"], linewidth=0.15, edgecolor="black", alpha=0.5)
    sns.kdeplot(ax=ax[2], data=stagnant_list, x="ti_trans", y="time_interval_trans", levels=10, color="brown", linewidths=linewidth, alpha=0.5, fill=True)
    sns.rugplot(ax=ax[2], data=stagnant_list, x="ti_trans", y="time_interval_trans", hue = "lithology", palette=colors_tfin, linewidth=rugline, alpha=rugalpha)
    # sns.regplot(ax=ax[2], data=stagnant_list, x="ti_trans", y="time_interval_trans", scatter=False, color="red", line_kws={"linewidth": linewidth})

    # ax2 = ax[2].twiny()
    # sns.histplot(ax=ax2, data=stagnant_list, y="time_interval_kin", hue = "lithology", palette=colors_tfin,  bins=20, alpha=0.3, kde=False, zorder = 1)
    # sns.histplot(ax=ax2, data=stagnant_list, y="time_interval_dyn", hue = "lithology", palette=colors_tfin,  bins=15, alpha=0.3, kde=False, zorder = 2)
    # sns.histplot(ax=ax2, data=stagnant_list, y="time_interval_trans", hue = "lithology", palette=colors_tfin,  bins=2, alpha=0.3, kde=False, zorder = 3)


    # Set plot limits and labels
    ax[2].set_xlim(0, 50)
    ax[2].set_ylim(0, 35)
    ax[2].axvline(x=35, color="grey", linestyle="--", linewidth=linewidth)
    ax[2].set_ylabel("Duration (Ma)")
    ax[2].set_xlabel("Time (Ma)")

    # Adjust layout and make last row shorter
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.)  # Adjust space between plots

    plt.savefig(f"{plot_loc}/time_vs_duration_scatter.png", dpi = 500)
    plt.close()

if __name__ == "__main__":
    main()