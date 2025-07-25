#!/usr/bin/env python3

import pandas as pd
import numpy as np
import json
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import os


def main():
    # Parse the arguments
    parser = argparse.ArgumentParser(description='Script that gets some models and gives the kinematic indicators')
    parser.add_argument('json_file', help='json file with model name, time at the end of computation, folder where to save the plots, legend')
    args = parser.parse_args()

    json_loc = '/home/vturino/PhD/projects/exhumation/pyInput/'
    json_file = f"{json_loc}{args.json_file}"

    # Read the json file
    with open(f"{json_file}") as json_file:
        configs = json.load(json_file)

    colors_tfin = {
        "sed": "cornflowerblue",
        "oc": "darkorange",
        "ecl": "forestgreen",
        "serp": "lightsalmon"
    }

    cr_colors = [
        "maroon",
        "red",
        "orange"
    ]

    alpha = 0.5
    linewidth = 1

    figloc = f"/home/vturino/PhD/projects/exhumation/plots/comparison/"
    plotname = f"Exhumation_threshold_velocity.pdf"
    count = len(configs["models"])

    time_max = 65

    # Create the plot with one column
    fig, ax = plt.subplots(count + 1, 1, figsize=(6, 7), height_ratios=[0.35] + [1] * count)

    for ind_m, m in enumerate(configs["models"]):
        # Create the folders to save the plots
        plot_loc = f"/home/vturino/PhD/projects/exhumation/plots/single_models/{m}"
        txt_loc = f'{plot_loc}/txt_files'
        if not os.path.exists(txt_loc):
            os.mkdir(txt_loc)

        # Load the data
        exhumed_list = pd.read_csv(f"{txt_loc}/exhumed_times.txt", sep="\s+")

        cr = pd.read_csv(f"{txt_loc}/2D_v.txt", sep="\s+")
        cr.iloc[0] = cr.iloc[1]
        cr["conv_rate"].iloc[0] = np.nan
        ax[0].plot(cr["time"] / 1e6, cr["conv_rate"], color=cr_colors[ind_m], linewidth=1)
        ax[0].axvline(x=35, color='grey', linestyle="--", linewidth=linewidth)
        ax[0].label_outer()
        ax[0].set_ylim(0, cr["conv_rate"].max() + 1)
        ax[0].set_xlim(13, 50)
        ax[0].set_title("Exhumable particles", fontsize=12)
        ax[0].set_yticks([0, 5])
        ax[0].set_xticks([])

        # Histograms for exhumed particles
        bin_edges = np.arange(0, time_max, 1)
        exhumed_list_sorted = exhumed_list.sort_values(by='lithology', ascending=True)
        sns.histplot(exhumed_list_sorted, x="texh", hue="lithology", bins=bin_edges, ax=ax[ind_m + 1],
                     palette=colors_tfin, alpha=alpha, linewidth=1, element="step", legend=False)

        # Set scale, labels, and limits
        ax[ind_m + 1].set_yscale("log")
        ax[ind_m + 1].label_outer()
        ax[ind_m + 1].set_xlabel("")
        ax[ind_m + 1].set_ylabel("Particles count (log)")
        ax[ind_m + 1].set_xlim(13, 50)
        ax[ind_m + 1].set_xticks([])
        ax[ind_m + 1].set_ylim(0.7, 5.e3)
        ax[ind_m + 1].axvline(x=35, color="grey", linestyle="--", label="35 Ma", linewidth=linewidth)

    ax[count].set_xticks([20, 30, 40, 50])
    ax[count].set_xlabel("Time (Myr)")

    plt.tight_layout()
    plt.subplots_adjust(hspace=0)

    plt.savefig(f"{figloc}/{plotname}", dpi=500)
    plt.close()


if __name__ == "__main__":
    main()
