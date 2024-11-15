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
    exhumed_list = pd.read_csv(f"{txt_loc}/timing_exhumed_particles.txt", sep="\s+")
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
    alpha = 1
    edgecolors_exh = "k"
    edgecolors_stag = "k"

    # #####################
    # ###### PLOT 1 #######
    # #####################

    # # Create the plot
    # fig, ax = plt.subplots(2, 1, figsize=(15, 10), gridspec_kw={'height_ratios': [1, 2]})

    # # Plot 0: convergence rate
    # ax[0].plot(cr["time"]/1e6, cr["conv_rate"], color="grey", label="Convergence rate", linewidth=3)
    # ax[0].set_ylabel("Convergence rate (cm/yr)", color="black", fontsize=labels_font)
    # ax[0].set_xlabel("Time (Ma)", fontsize=labels_font)
    # ax[0].set_title("Convergence rate and exhumed particles", fontsize=15)
    # ax[0].set_ylim(0, cr["conv_rate"].max()+0.2)
    # ax[0].set_xlim(0, 51)
    # ax[0].label_outer()  # Only show outer labels and tick labels
    # # vertical line at 35 Ma
    # ax[0].axvline(x=35, color="grey", linestyle="--", label="35 Ma")

    # # Plot 1: Pressure-time for exhumed particles
    # ax[1].add_patch(plt.Rectangle((0, 0), exhumed_list["tin"].max(), 2000, color='grey', alpha=alpha))
    
    # # Create legend handles manually for tfin and tmax using unique lithologies and colors
    # unique_lithologies = exhumed_list["lithology"].unique()
    # handles_tmax = [mpatches.Patch(color=colors_tmax[lith], label=lith) for lith in unique_lithologies]
    # handles_tfin = [mpatches.Patch(color=colors_tfin[lith], label=lith) for lith in unique_lithologies]



    # # Histograms for tmax and tfin (Manual layering with plt.bar)
    # bin_edges = np.arange(0, 50.5, 0.5)  # Define bin edges
    # width = np.diff(bin_edges)

    # # Plot tmax bars (lower count in front, higher in back)
    # for i in range(len(bin_edges) - 1):
    #     counts = [
    #         (np.histogram(exhumed_list[exhumed_list['lithology'] == lith]['tmax'], bins=bin_edges)[0], colors_tmax[lith]) 
    #         for lith in unique_lithologies
    #     ]
        
    #     # Sort counts by value so that the smaller count bars are drawn in front
    #     sorted_counts = sorted(counts, key=lambda x: x[0][i])
        
    #     # Plot the smaller count bar (in front)
    #     ax[1].bar(bin_edges[i], sorted_counts[0][0][i], width=width[i], color=sorted_counts[0][1],
    #               edgecolor=edgecolors_exh, linewidth=0.5, zorder=3)
    #     # Plot the larger count bar (in back)
    #     ax[1].bar(bin_edges[i], sorted_counts[1][0][i], width=width[i], color=sorted_counts[1][1],
    #               edgecolor=edgecolors_exh, linewidth=0.5, zorder=1)

    # # Plot tfin bars (lower count in front, higher in back)
    # for i in range(len(bin_edges) - 1):
    #     counts = [
    #         (np.histogram(exhumed_list[exhumed_list['lithology'] == lith]['tfin'], bins=bin_edges)[0], colors_tfin[lith]) 
    #         for lith in unique_lithologies
    #     ]
        
    #     # Sort counts by value so that the smaller count bars are drawn in front
    #     sorted_counts = sorted(counts, key=lambda x: x[0][i])
        
    #     # Plot the smaller count bar (in front)
    #     ax[1].bar(bin_edges[i], sorted_counts[0][0][i], width=width[i], color=sorted_counts[0][1],
    #               edgecolor=edgecolors_exh, linewidth=0.5, zorder=3)
    #     # Plot the larger count bar (in back)
    #     ax[1].bar(bin_edges[i], sorted_counts[1][0][i], width=width[i], color=sorted_counts[1][1],
    #               edgecolor=edgecolors_exh, linewidth=0.5, zorder=1)

    # # Add custom legends for tmax and tfin
    # legend_tmax = ax[1].legend(handles=handles_tmax, title="Peak Pressure", loc=(0.75, 0.85), frameon=True)
    # ax[1].add_artist(legend_tmax)  # Manually add the tmax legend
    # ax[1].legend(handles=handles_tfin, title="Time at Exhumability", loc=(0.85, 0.85), frameon=True)

    # # Set scale, labels, and limits
    # ax[1].set_yscale("log")
    # ax[1].text(1, 20, " Particles \n Subduction", fontsize=18)
    # ax[1].set_xlabel("Time (Ma)", fontsize=labels_font)
    # ax[1].set_ylabel("Particles count (log)", fontsize=labels_font)
    # ax[1].set_xlim(0, 51)
    # ax[1].set_ylim(1, 2.e3)
    # ax[1].axvline(x=35, color="grey", linestyle="--", label="35 Ma")

    # # Add no vertical space between the plots
    # plt.subplots_adjust(hspace=0)

    # plt.savefig(f"{plot_loc}/exhumVScr_hist.png", dpi=500)
    # plt.close()







    #####################
    ###### PLOT 2 #######
    #####################

    # ax[0]: convergence rate
    # ax[1]: Count of particles for each 0.5 Myr for exhumed particles, tfin
    # ax[2]: Count of particles for each 0.5 Myr for stagnant particles, tm_kin, tm_dyn, tm_trans

    # # Create the plot
    fig, ax = plt.subplots(3, 1, figsize=(17, 15), gridspec_kw={'height_ratios': [0.5, 1.5, 1.5]})

    # Plot 0: convergence rate
    ax[0].plot(cr["time"]/1e6, cr["conv_rate"], color="grey", label="Convergence rate (cm/yr)", linewidth=3)
    ax[0].axvline(x=35, color="grey", linestyle="--")
    ax[0].label_outer()  # Only show outer labels and tick labels
    ax[0].legend(loc='upper right')
    ax[0].set_ylim(0, cr["conv_rate"].max()+0.3)
    ax[0].set_xlim(0, 50)


    # Plot 1: Count of particles for each 0.5 Myr for exhumed particles, tfin, do as in first plot
    ax[1].add_patch(plt.Rectangle((0, 0), exhumed_list["tin"].max(), 2000, color='grey', alpha=alpha))
    
    # Create legend handles manually for tfin and tmax using unique lithologies and colors
    unique_lithologies_exhumed = exhumed_list["lithology"].unique()
    handles_tfin_exhumed = [mpatches.Patch(color=colors_tfin[lith], label=lith) for lith in unique_lithologies_exhumed]



    # Histograms for tmax and tfin (Manual layering with plt.bar)
    bin_edges = np.arange(0, 50.5, 0.5)  # Define bin edges
    width = np.diff(bin_edges)

    # Plot tfin bars (lower count in front, higher in back)
    for i in range(len(bin_edges) - 1):
        counts = [
            (np.histogram(exhumed_list[exhumed_list['lithology'] == lith]['tfin'], bins=bin_edges)[0], colors_tfin[lith]) 
            for lith in unique_lithologies_exhumed
        ]
        
        # Sort counts by value so that the smaller count bars are drawn in front
        sorted_counts = sorted(counts, key=lambda x: x[0][i])
        
        # Plot the smaller count bar (in front)
        ax[1].bar(bin_edges[i], sorted_counts[0][0][i], width=width[i], color=sorted_counts[0][1],
                  edgecolor=edgecolors_exh, linewidth=0.5, zorder=3, alpha=1)
        # Plot the larger count bar (in back)
        ax[1].bar(bin_edges[i], sorted_counts[1][0][i], width=width[i], color=sorted_counts[1][1],
                  edgecolor=edgecolors_exh, linewidth=0.5, zorder=1, alpha=1)

    # Add custom legends for tfin
    ax[1].legend(handles=handles_tfin_exhumed, title="Time at Exhumability", loc=(0.85, 0.85), frameon=True)

    # Set scale, labels, and limits
    ax[1].set_yscale("log")
    ax[1].label_outer()  # Only show outer labels and tick labels
    ax[1].text(1, 20, " Particles \n Subduction", fontsize=18)
    ax[1].set_xlabel("Time (Ma)", fontsize=labels_font)
    ax[1].set_ylabel("Particles count (log)", fontsize=labels_font)
    ax[1].set_xlim(0, 50)
    ax[1].set_ylim(0.9, 2.e3)
    ax[1].axvline(x=35, color="grey", linestyle="--", label="35 Ma")


    # Plot 2: Count of particles for each 0.5 Myr for stagnant particles, tm_kin, tm_dyn, tm_trans: do like for exhumed
    ax[2].add_patch(plt.Rectangle((0, 0), exhumed_list["tin"].max(), 2000, color='grey', alpha=alpha))

    unique_lithologies_stagnant = stagnant_list["lithology"].unique()
    handles_stagnant = [mpatches.Patch(color=colors_tfin[lith], label=lith) for lith in unique_lithologies_stagnant]

    # Define bin edges for tm_kin
    bin_edges = np.arange(0, 50, 0.5)  # Bin edges from 0 to 50 with step 0.5
    width = np.diff(bin_edges)

    # Initialize a dictionary to store the histogram counts for each lithology
    histograms = {}

    # Calculate histograms for tm_kin for each lithology
    for lith in unique_lithologies_stagnant:
        counts, _ = np.histogram(stagnant_list[stagnant_list['lithology'] == lith]['tm_kin'], bins=bin_edges)
        histograms[lith] = counts

    # Initialize a list to store (count, lithology) pairs
    sorted_histograms = []

    # Add all lithology counts to the sorted list
    for lith in unique_lithologies_stagnant:
        sorted_histograms.append((histograms[lith], lith))

    # Sort the histograms by counts (smaller counts will be first)
    sorted_histograms.sort(key=lambda x: np.sum(x[0]))  # Sort by the sum of counts for each lithology

    # Initialize an array for the cumulative counts (to stack the bars)
    cumulative_counts = np.zeros_like(bin_edges[:-1])

    # Plot the histograms, ensuring that smaller counts are plotted first
    for counts, lith in sorted_histograms:
        ax[2].bar(bin_edges[:-1], counts, width=width, label=lith, color=colors_tfin[lith], 
                alpha=1., bottom=cumulative_counts, edgecolor=edgecolors_stag, linewidth = 0.5)  # "bottom" stacks the bars
        cumulative_counts += counts  # Update cumulative counts


    # Now this but for tm_dyn
    # Initialize a dictionary to store the histogram counts for each lithology
    histograms = {}

    # Calculate histograms for tm_kin for each lithology
    for lith in unique_lithologies_stagnant:
        counts, _ = np.histogram(stagnant_list[stagnant_list['lithology'] == lith]['tm_dyn'], bins=bin_edges)
        histograms[lith] = counts

    # Initialize a list to store (count, lithology) pairs
    sorted_histograms = []

    # Add all lithology counts to the sorted list
    for lith in unique_lithologies_stagnant:
        sorted_histograms.append((histograms[lith], lith))

    # Sort the histograms by counts (smaller counts will be first)
    sorted_histograms.sort(key=lambda x: np.sum(x[0]))  # Sort by the sum of counts for each lithology

    # Initialize an array for the cumulative counts (to stack the bars)
    cumulative_counts = np.zeros_like(bin_edges[:-1])

    # Plot the histograms, ensuring that smaller counts are plotted first
    for counts, lith in sorted_histograms:
        ax[2].bar(bin_edges[:-1], counts, width=width, label=lith, color=colors_tfin[lith], 
                alpha=1., bottom=cumulative_counts, edgecolor=edgecolors_stag, linewidth = 0.5)  # "bottom" stacks the bars
        cumulative_counts += counts  # Update cumulative counts





    # Now this but for tm_trans
    # Initialize a dictionary to store the histogram counts for each lithology
    histograms = {}

    # Calculate histograms for tm_kin for each lithology
    for lith in unique_lithologies_stagnant:
        counts, _ = np.histogram(stagnant_list[stagnant_list['lithology'] == lith]['tm_trans'], bins=bin_edges)
        histograms[lith] = counts


    # Initialize a list to store (count, lithology) pairs
    sorted_histograms = []

    # Add all lithology counts to the sorted list
    for lith in unique_lithologies_stagnant:
        sorted_histograms.append((histograms[lith], lith))

    # Sort the histograms by counts (smaller counts will be first)
    sorted_histograms.sort(key=lambda x: np.sum(x[0]))  # Sort by the sum of counts for each lithology

    # Initialize an array for the cumulative counts (to stack the bars)
    cumulative_counts = np.zeros_like(bin_edges[:-1])

    # Plot the histograms, ensuring that smaller counts are plotted first
    for counts, lith in sorted_histograms:
        ax[2].bar(bin_edges[:-1], counts, width=width, label=lith, color=colors_tfin[lith], 
                alpha=1., bottom=cumulative_counts, edgecolor=edgecolors_stag, linewidth = 0.5)  # "bottom" stacks the bars
        cumulative_counts += counts  # Update cumulative counts









    # Add a title, labels, and legend
    ax[2].text(1, 20, " Particles \n Subduction", fontsize=18)
    ax[2].legend(handles=handles_stagnant, title="Lithologies")
    ax[2].set_xlim(0, 50)
    ax[2].set_ylim(0.9, 2.e3)
    ax[2].set_yscale("log")
    ax[2].axvline(x=35, color="grey", linestyle="--", label="35 Ma")



   

    # Add no vertical space between the plots
    plt.subplots_adjust(hspace=0)

    plt.savefig(f"{plot_loc}/CR_exhumed_stagnant.png", dpi=500)
    plt.close()





if __name__ == "__main__":
    main()
