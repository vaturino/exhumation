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


    alpha = 0.5
    rugline = 2
    rugalpha = 0.1
    linewidth = 1
    msize = 50

 
    # # Create the plot
    fig,ax = plt.subplots(3, 2, figsize=(15, 6), height_ratios=[0.25, 1, 1])

    # Plot 0: convergence rate
    ax[0,0].plot(cr["time"]/1e6, cr["conv_rate"], color="grey", label="Convergence rate (cm/yr)", linewidth=1)
    ax[0,0].axvline(x=35, color="grey", linestyle="--", linewidth = linewidth)
    ax[0,0].label_outer()  # Only show outer labels and tick labels
    ax[0,0].legend(loc='upper right')
    ax[0,0].set_ylim(0, cr["conv_rate"].max()+1)
    ax[0,0].set_xlim(0, 50)
    ax[0,0].set_title("Exhumable particles", fontsize=12)
    ax[1,0].set_yticks([0, 5])

    # Histograms for tmax and tfin (Manual layering with plt.bar)
    bin_edges = np.arange(0, 50,1)  # Define bin edges


    # Sort exhumed_list by lithology to ensure 'oc' is plotted over 'sed'
    exhumed_list_sorted = exhumed_list.sort_values(by='lithology', ascending=True)
    sns.histplot(exhumed_list_sorted, x="ti", hue="lithology", bins=bin_edges, ax=ax[1,0], palette=colors_tfin, alpha=alpha, linewidth=1, element="step", legend=False)
    # sns.histplot(exhumed_list_sorted, x="tm", hue="lithology", bins=bin_edges, ax=ax[1,0], palette=colors_tfin, alpha=alpha, linewidth=1, element="step", legend = False)
    

    # Set scale, labels, and limits
    ax[1,0].set_yscale("log")
    ax[1,0].label_outer()  # Only show outer labels and tick labels
    # ax[1,0].text(1, 20, " Particles \n Subduction", fontsize=18)
    ax[1,0].set_xlabel("")
    ax[1,0].set_ylabel("Particles count (log)")
    ax[1,0].set_xlim(0, 50)
    ax[1,0].set_ylim(0.7, 5.e3)
    ax[1,0].axvline(x=35, color="grey", linestyle="--", label="35 Ma", linewidth = linewidth)
    ax[1,0].text(0.32, 0.9, "Time at exhumable threshold", fontsize=12, transform=ax[1,0].transAxes)


    # Add pathch forom the first bin (smaller x value) to end of the plot

        # Plot Exhumed data
    sns.scatterplot(ax=ax[2,0], data=exhumed_list, x="ti", y="time_interval", hue="lithology", palette=colors_tfin, zorder=10, 
                    style="lithology", s = msize, markers = ["o", "d", "s", "X"], linewidth=0.15, edgecolor="black", alpha=0.8, legend=False)
    # # Rugplot with variable height using ax.vlines
    # for _, row in exhumed_list.iterrows():
    #     ax[2,0].hlines(
    #         y = row["time_interval"],
    #         xmin = row["ti"],
    #         xmax = row["ti"] + row["time_interval"],
    #         color = colors_tfin[row["lithology"]],
    #         linewidth = rugline,
    #         alpha = rugalpha
    #     )
    sns.rugplot(ax=ax[2,0], data=exhumed_list, x=None, y="time_interval", hue = "lithology", palette=colors_tfin, linewidth=rugline, alpha=rugalpha, height=0.1)
    ax[2,0].set_xlim(0, 50)
    ax[2,0].axvline(x=35, color="grey", linestyle="--", linewidth=linewidth)
    ax[2,0].set_ylabel("Duration (Myr)")
    ax[2,0].set_xlabel("Time (Myr)")
    ax[2,0].set_ylim(0, 35)
    ax[2,0].set_yticks([0, 5, 10, 15, 20, 25, 30, 35])


    # COLUMN 2

    # Plot 0: convergence rate
    ax[0,1].plot(cr["time"]/1e6, cr["conv_rate"], color="grey", label="Convergence rate (cm/yr)", linewidth=1)
    ax[0,1].axvline(x=35, color="grey", linestyle="--", linewidth = linewidth)
    ax[0,1].label_outer()  # Only show outer labels and tick labels
    ax[0,1].legend(loc='upper right')
    ax[0,1].set_ylim(0, cr["conv_rate"].max()+1)
    ax[0,1].set_xlim(0, 50)
    ax[0,1].set_title("Stagnant particles", fontsize=12)
    ax[0,1].set_yticks([0, 5])



    # Histograms for tmax and tfin (Manual layering with plt.bar)
    bin_edges = np.arange(0, 50, 1)  # Define bin edges
    # Histograms for tmax and tfin (Manual layering with plt.bar)
    stagnant_list["ti"] = np.nan

    # Handle ti_kin, ti_dyn, ti_trans
    new_rows = []
    for index, row in stagnant_list.iterrows():
        ti_values = [row["ti_kin"], row["ti_dyn"], row["ti_trans"]]
        non_nan_ti_values = [ti for ti in ti_values if not np.isnan(ti)]
        if len(non_nan_ti_values) == 1:
            stagnant_list.at[index, "ti"] = non_nan_ti_values[0]
        elif len(non_nan_ti_values) > 1:
            for ti in non_nan_ti_values:
                new_row = row.copy()
                new_row["ti"] = ti
                new_rows.append(new_row)
            stagnant_list.at[index, "ti"] = np.nan  # Mark original row as NaN to be dropped later

    # Concatenate new rows to the DataFrame and drop original rows with NaN "ti"
    stagnant_list = pd.concat([stagnant_list, pd.DataFrame(new_rows)], ignore_index=True)
    stagnant_list = stagnant_list.dropna(subset=["ti"])

    sns.histplot(stagnant_list, x="ti", ax=ax[1,1], hue = "lithology", palette= colors_tfin, alpha=alpha, linewidth=1, element="step", bins= bin_edges)
    # stagnant_list_sorted = stagnant_list.sort_values(by='lithology', ascending=True)

    # sns.histplot(stagnant_list_sorted, x="ti_kin", hue="lithology", ax=ax[1,1], palette=colors_tfin, alpha=alpha, linewidth=1, element="step", bins= bin_edges, legend = False)
    # sns.histplot(stagnant_list_sorted, x="ti_dyn", hue="lithology", ax=ax[1,1], palette=colors_tfin, alpha=alpha, linewidth=1, element="step", bins= bin_edges, legend = False)
    # sns.histplot(stagnant_list_sorted, x="ti_trans", hue="lithology", ax=ax[1,1], palette=colors_tfin, alpha=alpha, linewidth=1, element="step", bins= bin_edges, legend=True)

    ax[1,1].text(17, 2500, "Beginning of stagnation", fontsize=12)

    # Add a title, labels, and legend
    ax[1,1].set_xlabel("")
    ax[1,1].set_ylabel("")
    ax[1,1].set_xlim(0, 50)
    ax[1,1].set_ylim(0.7, 5.e3)
    ax[1,1].set_xticks([])
    #put y ticks to the right
    # ax[1,1].yaxis.tick_right()
    ax[1,1].set_yscale("log")
    ax[1,1].axvline(x=35, color="grey", linestyle="--", label="35 Ma", linewidth = linewidth)


    # Plot initial time and duration of stagnant particles
    sns.scatterplot(ax=ax[2,1], data=stagnant_list, x="ti_dyn", y="time_interval_dyn", hue="lithology", palette=colors_tfin, zorder=10, 
                    style="lithology", s = msize, markers = ["o", "d", "s", "X"], linewidth=0.15, edgecolor="black", alpha=0.8)
    sns.scatterplot(ax=ax[2,1], data=stagnant_list, x="ti_kin", y="time_interval_kin", hue="lithology", palette=colors_tfin, zorder=10, 
                    style="lithology", s = msize, markers = ["o", "d", "s", "X"], linewidth=0.15, edgecolor="black", alpha=0.8, legend=False)
    sns.scatterplot(ax=ax[2,1], data=stagnant_list, x="ti_trans", y="time_interval_trans", hue="lithology", palette=colors_tfin, zorder=10, 
                    style="lithology", s = msize, markers = ["o", "d", "s", "X"], linewidth=0.15, edgecolor="black", alpha=0.8, legend=False)  
    ax[2,1].legend(loc='upper left') 
 

    sns.rugplot(ax=ax[2,1], data=stagnant_list, x=None, y="time_interval_dyn", hue = "lithology", palette=colors_tfin, linewidth=rugline, alpha=rugalpha, height=0.1)
    sns.rugplot(ax=ax[2,1], data=stagnant_list, x=None, y="time_interval_kin", hue = "lithology", palette=colors_tfin, linewidth=rugline, alpha=rugalpha, height=0.1)
    sns.rugplot(ax=ax[2,1], data=stagnant_list, x=None, y="time_interval_trans", hue="lithology", palette=colors_tfin, linewidth=rugline, alpha=rugalpha, height=0.1)
    
    # for _, row in stagnant_list.iterrows():
    #     ax[2,1].hlines(
    #         y = row["time_interval_dyn"],
    #         xmin = row["ti_dyn"],
    #         xmax = row["ti_dyn"] + row["time_interval_dyn"],
    #         color = colors_tfin[row["lithology"]],
    #         linewidth = rugline,
    #         alpha = rugalpha
    #     )
    
    # for _, row in stagnant_list.iterrows():
    #     ax[2,1].hlines(
    #         y = row["time_interval_kin"],
    #         xmin = row["ti_kin"],
    #         xmax = row["ti_kin"] + row["time_interval_kin"],
    #         color = colors_tfin[row["lithology"]],
    #         linewidth = rugline,
    #         alpha = rugalpha
    #     )

    # for _, row in stagnant_list.iterrows():
    #     ax[2,1].hlines(
    #         y = row["time_interval_trans"],
    #         xmin = row["ti_trans"],
    #         xmax = row["ti_trans"] + row["time_interval_trans"],
    #         color = colors_tfin[row["lithology"]],
    #         linewidth = rugline,
    #         alpha = rugalpha
    #     )


    ax[2,1].set_xlim(0, 50)
    ax[2,1].axvline(x=35, color="grey", linestyle="--", linewidth=linewidth)
    ax[2,1].set_ylabel("")
    ax[2,1].set_xlabel("Time (Myr)")
    ax[2,1].set_ylim(0, 35)
    # ax[2,1].yaxis.tick_right()
    # ax[2,1].set_yticks([0, 5, 10, 15, 20, 25, 30, 35])



   

    # Add no vertical space between the plots
    plt.tight_layout()
    plt.subplots_adjust(hspace=0, wspace=0.06)

    plt.savefig(f"{plot_loc}/CR_timing_exhumed_rugplot.png", dpi=500)
    plt.close()





if __name__ == "__main__":
    main()
