#!/usr/bin/env python3

import pandas as pd
import numpy as np
import json
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import os
import matplotlib.patches as mpatches
from scipy.stats import gaussian_kde


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

    # Create the folders to save the plots
    plot_loc = f"/home/vturino/PhD/projects/exhumation/plots/single_models/{configs['models'][0]}"
    txt_loc = f'{plot_loc}/txt_files'
    if not os.path.exists(txt_loc):
        os.mkdir(txt_loc)

    # Load the data
    exhumed_list = pd.read_csv(f"{txt_loc}/exhumed_times.txt", sep="\s+")

    

    cr = pd.read_csv(f"{txt_loc}/2D_v.txt", sep="\s+")
    cr["conv_rate"].iloc[0] = np.nan

    xmax_plot = 50
    time_max = 65

    colors_tfin = {
        "sed": "cornflowerblue",
        "oc": "darkorange",
        "ecl": "#61aa2f",
        "serp": "lightsalmon"
    }

    colors_tin = {
        "sed": "mediumblue",
        "oc": "chocolate",
        "ecl": "darkgreen",
        "serp": "orangered"
    }


    alpha = 0.5
    linewidth = 1

 
    # # Create the plot
    fig,ax = plt.subplots(5, 1, figsize=(10, 10), height_ratios=[0.25, 0.05, 1, 0.05, 1])

    # Plot 0: convergence rate
    ax[0].plot(cr["time"]/1e6, cr["conv_rate"], color="grey", label="Convergence rate (cm/yr)", linewidth=1)
    ax[0].axvline(x=35, color="grey", linestyle="--", linewidth = linewidth)
    ax[0].label_outer()  # Only show outer labels and tick labels
    ax[0].legend(loc='upper right')
    ax[0].set_ylim(0, cr["conv_rate"].max()+1)
    ax[0].set_xlim(0, xmax_plot)
    ax[0].set_title("Exhumable particles", fontsize=12)
    ax[0].set_yticks([0, 5])

    ax[1].set_visible(False)

    # # Histograms for tmax and tfin (Manual layering with plt.bar)
    bin_edges = np.arange(0, time_max,1)  # Define bin edges


    # Sort exhumed_list by lithology to ensure 'oc' is plotted over 'sed'
    exhumed_list_sorted = exhumed_list.sort_values(by='lithology', ascending=True)
    sns.histplot(exhumed_list_sorted, x="tin", hue="lithology", bins=bin_edges, ax=ax[2], palette=colors_tfin, alpha=alpha, linewidth=1, element="step", legend=False)
    

    # Set scale, labels, and limits
    ax[2].set_yscale("log")
    ax[2].label_outer()  # Only show outer labels and tick labels
    # ax[2].text(1, 20, " Particles \n Subduction", fontsize=18)
    ax[2].set_xlabel("")
    ax[2].set_ylabel("Particles count (log)")
    ax[2].set_xlim(0, xmax_plot)
    ax[2].set_ylim(0.7, 5.e3)
    ax[2].axvline(x=35, color="grey", linestyle="--", label="35 Ma", linewidth = linewidth)
    ax[2].text(0.32, -1, "Exhumation onset", fontsize=12, transform=ax[1].transAxes)


    ax[3].set_visible(False)


    sns.histplot(exhumed_list_sorted, x="texh", hue="lithology", bins=bin_edges, ax=ax[4], palette=colors_tfin, alpha=0.5, linewidth=1, element="step", legend=False)
    ax[4].set_yscale("log")
    ax[4].set_ylabel("Particles count (log)")
    ax[4].set_xlabel("Time (Myr)")
    ax[4].set_xlim(0, xmax_plot)
    ax[4].set_ylim(0.7, 5.e3)
    ax[4].axvline(x=35, color="grey", linestyle="--", label="35 Ma", linewidth = linewidth)
    ax[4].text(0.32, -1, "Time at exhumable threshold", fontsize=12, transform=ax[3].transAxes)



   
    plt.subplots_adjust(hspace=0.0)





    

    plt.savefig(f"{plot_loc}/time_at_end.pdf", dpi=500)
    plt.close()






if __name__ == "__main__":
    main()
