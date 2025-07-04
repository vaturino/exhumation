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

    plotloc = f'/home/vturino/PhD/projects/exhumation/plots/{configs["plot_folder"][0]}'
    plotname = "rheology_comparison.pdf"

    crfolder = f"/home/vturino/PhD/projects/exhumation/plots/single_models/{configs['models'][0]}/txt_files"

    cr = pd.read_csv(f"{crfolder}/2D_v.txt", sep="\s+")
    cr["conv_rate"].iloc[0] = np.nan


    # colors
    color = ["mediumslateblue", "darkslateblue", "forestgreen", "darkgreen"]



    fig, ax = plt.subplots(6, 2, figsize=(17, 10), height_ratios=[0.35, 0.3,  1, 1, 1, 1])
    
    # Plot 0: convergence rate
    ax[0,0].plot(cr["time"]/1e6, cr["conv_rate"], color="grey", label="Convergence rate (cm/yr)", linewidth=1)
    ax[0,0].axvline(x=35, color="grey", linestyle="--", linewidth = 1)
    # ax[0,0].label_outer()  # Only show outer labels and tick labels
    ax[0,0].set_xlim(0, 50)
    ax[0,0].set_ylim(0, cr["conv_rate"].max()+1)

    

    ax[0,1].plot(cr["time"]/1e6, cr["conv_rate"], color="grey", label="Convergence rate (cm/yr)", linewidth=1)
    ax[0,1].axvline(x=35, color="grey", linestyle="--", linewidth = 1)
    # ax[0,1].label_outer()  # Only show outer labels and tick labels
    ax[0,1].set_xlim(0, 50)
    ax[0,1].set_ylim(0, cr["conv_rate"].max()+1)

    ax[1,0].set_visible(False)
    ax[1,1].set_visible(False)




    for m, mod in enumerate(configs['models']):
        mod_loc = f"/home/vturino/PhD/projects/exhumation/plots/single_models/{mod}"
        exhumed_list = pd.read_csv(f"{mod_loc}/txt_files/exhumed_times.txt", sep="\s+")
        stagnant_list = pd.read_csv(f"{mod_loc}/txt_files/stagnant_particles.txt", sep="\s+")

        bin_edges = np.arange(0, 50, 1)  # Define bin edges

        # Left hand side: exhumed particles
        ax[2, 0].set_title("Exhumed particles", fontsize=14)
        sns.histplot(
            exhumed_list,
            x="ti",
            ax=ax[m+2, 0],
            bins=bin_edges,
            color=color[m],
            alpha=0.5,
            linewidth=1,
            element="step",
            legend=True,
        )
        ax[m+2, 0].axvline(x=35, color="grey", linestyle="--", linewidth = 1)
        ax[m+2, 0].set_yscale("log")
        ax[m+2, 0].set_xlim(0, 50)
        ax[m+2, 0].set_ylim(0.7, 5.e3)
        ax[m+2, 0].set_ylabel("Count (log)")
        ax[m+2, 0].set_xticks([])
        ax[m+2, 0].set_xlabel("")
        ax[m+2, 0].text(0.1, 0.8, configs["names"][m], horizontalalignment='left', verticalalignment='center', transform=ax[m+2, 0].transAxes)


        # Right hand side: stagnant particles
        ax[2, 1].set_title("Stagnant particles", fontsize=14)

        # Handle ti_kin, ti_dyn, ti_trans
        new_rows = []
        for index, row in stagnant_list.iterrows():
            ti_values = [row["ti_dyn"], row["ti_kin"], row["ti_trans"]]
            time_interval_values = [row["time_interval_dyn"], row["time_interval_kin"], row["time_interval_trans"]]
            lithology_values = [row["lithology_dyn"], row["lithology_kin"], row["lithology_trans"]]
            
            for ti, time_interval, lithology in zip(ti_values, time_interval_values, lithology_values):
                if not np.isnan(ti):  # Only include non-NaN values
                    new_row = row.copy()
                    new_row["ti"] = ti
                    new_row["time_interval"] = time_interval
                    new_row["lith_time"] = lithology
                    new_rows.append(new_row)

        # Concatenate new rows to create the final DataFrame
        stagnant_list_expanded = pd.DataFrame(new_rows)

        
        sns.histplot(
            stagnant_list_expanded,
            x="ti",
            ax=ax[m+2, 1],
            bins=bin_edges,
            color=color[m],
            alpha=0.5,
            linewidth=1,
            element="step",
            legend=True,
        )

        ax[m+2, 1].axvline(x=35, color="grey", linestyle="--", linewidth = 1)
        ax[m+2, 1].set_yscale("log")
        ax[m+2, 1].set_xlim(0, 50)
        ax[m+2, 1].set_ylim(0.7, 5.e3)
        ax[m+2, 1].set_ylabel("Count (log)")
        ax[m+2, 1].set_xlabel("Time (Ma)")
        ax[m+2, 1].set_xticks([])
        ax[m+2, 1].set_xticklabels([])
        ax[m+2, 1].text(0.1, 0.8, configs["names"][m], horizontalalignment='left', verticalalignment='center', transform=ax[m+2, 1].transAxes)

    # Add x axis numbers just to the bottom row
    for i in range(2):
        ax[5, i].set_xlabel("Time (Ma)")
        ax[5, i].set_xticks(np.arange(0, 51, 10))
        ax[5, i].set_xticklabels(np.arange(0, 51, 10))




    
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.0, wspace=0.2)
    plt.savefig(f"{plotloc}/rheo_no_lith.pdf")
    plt.close()





if __name__ == "__main__":
    main()
