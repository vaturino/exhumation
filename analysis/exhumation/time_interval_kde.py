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
    width_exh = 0.3
    width_stag = 0.3
    linewidth = 1

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
    # hist_exh = sns.histplot(ax=ax[1], data=exhumed_list, x="ti", y="time_interval", bins=100, cbar=False,
                            # binwidth=(width_exh, 3*width_exh), color='magenta', cbar_kws={"label": "Exhumed", "shrink": 0.8}, zorder = 10)
    hist_exh = sns.histplot(ax=ax[1], data=exhumed_list, x="ti", y="time_interval", bins=100, cbar=False,
                            binwidth=(width_exh, 3*width_exh), hue = "lithology", zorder = 10, palette=colors_tfin)
    sns.kdeplot(ax=ax[1], data=exhumed_list, x="ti", y="time_interval", levels=5, linewidths=1, alpha=0.3, color="magenta", fill=True, label="Exhumed")
    ax[1].set_xlim(0, 50)
    ax[1].axvline(x=35, color="grey", linestyle="--", linewidth=linewidth)
    ax[1].set_xticks([])
    ax[1].set_ylabel("Duration (Ma)")
    ax[1].set_xlabel("")
    ax[1].set_ylim(0, 35)
    ax[1].set_yticks([5, 10, 15, 20, 25, 30, 35])
    # cbar_exh = fig.add_axes([0.05, 0.5, 0.015, 0.25])
    # cbar_exh = plt.colorbar(hist_exh.collections[0], cax=cbar_exh)

    # Plot Stagnant data (Dynamic slowdown)
    # hist_stag_dyn = sns.histplot(ax=ax[2], data=stagnant_list, x="ti_dyn", y="time_interval_dyn", bins=100, cbar=False,
    #                             binwidth=(width_stag, 3*width_stag), color='red', cbar_kws={"label": "Dynamic slowdown", "shrink": 0.8}, zorder = 10)
    hist_stag_dyn = sns.histplot(ax=ax[2], data=stagnant_list, x="ti_dyn", y="time_interval_dyn", bins=100, cbar=False,
                                binwidth=(width_stag, 3*width_stag), hue="lithology", zorder=10, palette=colors_tfin)
    # sns.kdeplot(ax=ax[2], data=stagnant_list, x="ti_dyn", y="time_interval_dyn", levels=5, linewidths=1, alpha=0.3, color="red", fill=True, label="Dynamic slowdown")


    # Plot Stagnant data (Kinematic slowdown)
    # hist_stag_kin = sns.histplot(ax=ax[2], data=stagnant_list, x="ti_kin", y="time_interval_kin", bins=100, cbar=False,
    #                             binwidth=(width_stag, 3*width_stag), color='green', cbar_kws={"label": "Kinematic slowdown", "shrink": 0.8})
    hist_stag_kin = sns.histplot(ax=ax[2], data=stagnant_list, x="ti_kin", y="time_interval_kin", bins=100, cbar=False,
                                binwidth=(width_stag, 3*width_stag), hue="lithology", zorder=10, palette=colors_tfin)
    sns.kdeplot(ax=ax[2], data=stagnant_list, x="ti_kin", y="time_interval_kin", levels=5, linewidths=1, alpha=0.3, color="green", fill=True, label="Kinematic slowdown")

    # Plot Stagnant data (Transition point)
    # hist_stag_trans = sns.histplot(ax=ax[2], data=stagnant_list, x="ti_trans", y="time_interval_trans", bins=100, cbar=False,
    #                             binwidth=(width_stag, 3*width_stag), color='blue', cbar_kws={"label": "Transition point", "shrink": 0.8})
    hist_stag_trans = sns.histplot(ax=ax[2], data=stagnant_list, x="ti_trans", y="time_interval_trans", bins=100, cbar=False,
                                binwidth=(width_stag, 3*width_stag), hue="lithology", zorder=10, palette=colors_tfin)
    sns.kdeplot(ax=ax[2], data=stagnant_list, x="ti_trans", y="time_interval_trans", levels=5, linewidths=1, alpha=0.3, color="blue", fill=True, label="Transition point")

    # cbar_dyn = fig.add_axes([0.05, 0.1, 0.015, 0.25])
    # cbar_dyn = plt.colorbar(hist_stag_dyn.collections[0], cax=cbar_dyn)

    # cbar_kin = fig.add_axes([0.15, 0.07, 0.01, 0.25])
    # cbar_kin = plt.colorbar(hist_stag_kin.collections[1], cax=cbar_kin)


    # Set plot limits and labels
    ax[2].set_xlim(0, 50)
    ax[2].set_ylim(0, 35)
    ax[2].axvline(x=35, color="grey", linestyle="--", linewidth=linewidth)
    ax[2].set_ylabel("Duration (Ma)")
    ax[2].set_xlabel("Time (Ma)")

    # Adjust layout and make last row shorter
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.)  # Adjust space between plots

    plt.savefig(f"{plot_loc}/time_vs_duration.png", dpi = 500)
    plt.close()

if __name__ == "__main__":
    main()