#!/usr/bin/env python3

import pandas as pd
import numpy as np
import json
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import os
import matplotlib.patches as mpatches


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
    stagnant_list = pd.read_csv(f"{txt_loc}/stagnant_particles.txt", sep="\s+")
    cr = pd.read_csv(f"{txt_loc}/2D_v.txt", sep="\s+")
    cr["conv_rate"].iloc[0] = np.nan



    # Create a mapping for lithology colors
     # Define color palettes
    colors_tin = {
        "sed": "midnightblue",
        "oc": "#733A11",
        "ecl": "#003300",
        "serp": "#3b0000"
    }

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


    labels_font = 13
    alpha = 0.5
    edgecolors_exh = "k"
    edgecolors_stag = "k"

   


    #####################
    ###### PLOT 2 #######
    #####################

    # ax[0]: convergence rate
    # ax[1]: Count of particles for each 0.5 Myr for exhumed particles, tfin
    # ax[2]: Count of particles for each 0.5 Myr for stagnant particles, tm_kin, tm_dyn, tm_trans

    alpha = 0.5

    # # Create the plot
    fig, ax = plt.subplots(3, 1, figsize=(8,5), gridspec_kw={'height_ratios': [0.5, 1.5, 1.5]})

    # Plot 0: convergence rate
    ax[0].plot(cr["time"]/1e6, cr["conv_rate"], color="grey", label="Convergence rate (cm/yr)", linewidth=1)
    ax[0].axvline(x=35, color="grey", linestyle="--", linewidth = 0.5)
    ax[0].label_outer()  # Only show outer labels and tick labels
    ax[0].legend(loc='upper right')
    ax[0].set_ylim(0, cr["conv_rate"].max()+1)
    ax[0].set_xlim(0, 50)

    # Histograms for tmax and tfin (Manual layering with plt.bar)
    bin_edges = np.arange(0, 50,1)  # Define bin edges


    # Sort exhumed_list by lithology to ensure 'oc' is plotted over 'sed'
    exhumed_list_sorted = exhumed_list.sort_values(by='lithology', ascending=True)
    sns.histplot(exhumed_list_sorted, x="ti", hue="lithology", bins=bin_edges, ax=ax[1], palette=colors_tfin, alpha=alpha, linewidth=1, element="step", legend=True)
    # sns.histplot(exhumed_list_sorted, x="tm", hue="lithology", bins=bin_edges, ax=ax[1], palette=colors_tfin, alpha=alpha, linewidth=1, element="step", legend = False)
    

    # Set scale, labels, and limits
    ax[1].set_yscale("log")
    ax[1].label_outer()  # Only show outer labels and tick labels
    # ax[1].text(1, 20, " Particles \n Subduction", fontsize=18)
    ax[1].set_xlabel("Time (Ma)")
    ax[1].set_ylabel("Particles count (log)")
    ax[1].set_xlim(0, 50)
    ax[1].set_ylim(0.9, 5.e3)
    ax[1].axvline(x=35, color="grey", linestyle="--", label="35 Ma", linewidth = 0.5)
    # Add pathch forom the first bin (smaller x value) to end of the plot
    




    # Histograms for tmax and tfin (Manual layering with plt.bar)
    bin_edges = np.arange(0, 50, 1)  # Define bin edges
    stagnant_list_sorted = stagnant_list.sort_values(by='lithology', ascending=True)

    sns.histplot(stagnant_list_sorted, x="ti_kin", hue="lithology", ax=ax[2], palette=colors_tfin, alpha=alpha, linewidth=1, element="step", bins= bin_edges, legend = False)
    sns.histplot(stagnant_list_sorted, x="ti_dyn", hue="lithology", ax=ax[2], palette=colors_tfin, alpha=alpha, linewidth=1, element="step", bins= bin_edges, legend = False)
    sns.histplot(stagnant_list_sorted, x="ti_trans", hue="lithology", ax=ax[2], palette=colors_tfin, alpha=alpha, linewidth=1, element="step", bins= bin_edges,    )




    # Add a title, labels, and legend
    # ax[2].text(1, 20, " Particles \n Subduction", fontsize=18)
    # ax[2].legend(handles=handles_stagnant, title="Lithologies")
    ax[2].set_xlim(0, 50)
    ax[2].set_ylim(0.9, 5.e3)
    ax[2].set_yscale("log")
    ax[2].axvline(x=35, color="grey", linestyle="--", label="35 Ma", linewidth = 0.5)



   

    # Add no vertical space between the plots
    plt.subplots_adjust(hspace=0)

    plt.savefig(f"{plot_loc}/CR_exhumed_stagnant.pdf", dpi=500)
    plt.close()





if __name__ == "__main__":
    main()
