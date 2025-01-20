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


    #colors for models
    colors = ["red", "grey", "blue"]

    alpha = 0.5
    rugline = 2
    rugalpha = 0.1
    linewidth = 1
    msize = 50

    figloc = f"/home/vturino/PhD/projects/exhumation/plots/rheology_comparison/{configs['test']}"
    # count number of models in configs["models"]
    n_models = len(configs["models"])

    # # Create the plot
    fig,ax = plt.subplots(n_models+1, 3, figsize=(20, n_models*3), height_ratios=[0.25] + [1] * n_models, width_ratios=[1, 1, 0.1])

    for ind_m, m in enumerate(configs["models"]):
        # Create the folders to save the plots
        plot_loc = f"/home/vturino/PhD/projects/exhumation/plots/single_models/{m}"
        txt_loc = f'{plot_loc}/txt_files'
        if not os.path.exists(txt_loc):
            os.mkdir(txt_loc)

        # Load the data
        exhumed_list = pd.read_csv(f"{txt_loc}/exhumed_times.txt", sep="\s+")
        stagnant_list = pd.read_csv(f"{txt_loc}/stagnant_particles.txt", sep="\s+")
        cr = pd.read_csv(f"{txt_loc}/2D_v.txt", sep="\s+")
        cr["conv_rate"].iloc[0] = np.nan

        # add stacked bar with amount of subducted, exhumed, and stagnant particles
        # Create a dictionary with stagnant, exhumed, and subducted particles
        classification = {
            "subducted": "lightcyan",
            "exhumed": "orchid",
            "stagnant": "yellowgreen"
        }

        #count amount of total particles:
        particles = pd.read_csv(f"{txt_loc}/particles_indexes.csv", sep="\s+")
        nparts = len(particles)
        #count amount of exhumed particles
        nexhumed = len(exhumed_list)
        #count amount of stagnant particles
        nstagnant = len(stagnant_list)

        #calculate percentage of particles exhumed/stagnant with respect to total
        pexhumed = nexhumed/nparts
        pstagnant = nstagnant/nparts
        psubducted = 1 - pexhumed - pstagnant


        # Plot 0: convergence rate
        if ind_m == 1:
            ax[0,0].plot(cr["time"]/1e6, cr["conv_rate"], color='grey', label="Convergence rate (cm/yr)", linewidth=1)
            ax[0,0].axvline(x=35, color="grey", linestyle="--", linewidth = linewidth)
            ax[0,0].label_outer()  # Only show outer labels and tick labels
            # ax[0,0].legend(loc='upper right')
            ax[0,0].set_ylim(0, cr["conv_rate"].max()+1)
            ax[0,0].set_xlim(0, 50)
            ax[0,0].set_title("Exhumable particles", fontsize=12)
            ax[0,0].set_yticks([0, 5])

        sns.scatterplot(ax=ax[ind_m+1,0], data=exhumed_list, x="ti", y="time_interval",  color=colors[ind_m], zorder=10, 
                        s = msize, markers = ["o", "d", "s", "X"], linewidth=0.15, edgecolor="black", alpha=0.5, legend=True, label=configs["legend"][ind_m])
        
        sns.rugplot(ax=ax[ind_m+1,0], data=exhumed_list, x=None, y="time_interval", color=colors[ind_m], zorder=10, linewidth=rugline, alpha=rugalpha)
        ax[ind_m+1,0].set_xlim(0, 50)
        ax[ind_m+1,0].axvline(x=35, color="grey", linestyle="--", linewidth=linewidth)
        ax[ind_m+1,0].set_ylim(0, 35)
        ax[ind_m+1,0].set_yticks([0, 10, 20, 20, 30])
        ax[ind_m+1,0].set_xlabel("")
        ax[ind_m+1,0].set_xticks([])
        ax[ind_m+1,0].set_ylabel("Duration (Ma)")


        #Col 2: Stagnant particles
        if ind_m == 1:
            ax[0,1].set_title("Stagnant particles", fontsize=12)
            ax[0,1].plot(cr["time"]/1e6, cr["conv_rate"], color="grey", label="Convergence rate (cm/yr)", linewidth=1)
            ax[0,1].axvline(x=35, color="grey", linestyle="--", linewidth = linewidth)
            ax[0,1].label_outer()
            ax[0,1].set_ylim(0, cr["conv_rate"].max()+1)
            ax[0,1].set_xlim(0, 50)

        sns.scatterplot(ax=ax[ind_m+1,1], data=stagnant_list, x="ti_dyn", y="time_interval_dyn", color=colors[ind_m], zorder=10, 
                        s = msize, markers = ["o", "d", "s", "X"], linewidth=0.15, edgecolor="black", alpha=0.5, legend=True)
        sns.scatterplot(ax=ax[ind_m+1,1], data=stagnant_list, x="ti_kin", y="time_interval_kin", color=colors[ind_m], zorder=10,
                        s = msize, markers = ["o", "d", "s", "X"], linewidth=0.15, edgecolor="black", alpha=0.5, legend=False)
        sns.scatterplot(ax=ax[ind_m+1,1], data=stagnant_list, x="ti_trans", y="time_interval_trans", color=colors[ind_m], zorder=10,
                        s = msize, markers = ["o", "d", "s", "X"], linewidth=0.15, edgecolor="black", alpha=0.5, legend=False)
        
        sns.rugplot(ax=ax[ind_m+1,1], data=stagnant_list, x=None, y="time_interval_dyn", color=colors[ind_m], zorder=10, linewidth=rugline, alpha=rugalpha)
        sns.rugplot(ax=ax[ind_m+1,1], data=stagnant_list, x=None, y="time_interval_kin", color=colors[ind_m], zorder=10, linewidth=rugline, alpha=rugalpha)
        sns.rugplot(ax=ax[ind_m+1,1], data=stagnant_list, x=None, y="time_interval_trans", color=colors[ind_m], zorder=10, linewidth=rugline, alpha=rugalpha)
        
        ax[ind_m+1,1].set_xlim(0, 50)
        ax[ind_m+1,1].axvline(x=35, color="grey", linestyle="--", linewidth=linewidth)
        ax[ind_m+1,1].set_ylabel("")
        ax[ind_m+1,1].set_xlabel("")
        ax[ind_m+1,1].set_ylim(0, 35)
        ax[ind_m+1,1].set_yticks([0, 10, 20, 20, 30])
        ax[ind_m+1,1].set_xticks([])
        ax[ind_m+1,1].set_ylabel("Duration (Ma)")


        #Col 3: vertically stacked bar charts, no axes
        bar_width = 0.05
        ax[ind_m+1,2].bar(0.8, psubducted, color=classification["subducted"], 
                        label="Subducted", edgecolor="black", linewidth=0.5, width=bar_width)
        ax[ind_m+1,2].bar(0.8, pexhumed, bottom=psubducted, color=classification["exhumed"], 
                        label="Exhumed", edgecolor="black", linewidth=0.5, width=bar_width)
        ax[ind_m+1,2].bar(0.8, pstagnant, bottom=psubducted+pexhumed, color=classification["stagnant"], 
                        label="Stagnant", edgecolor="black", linewidth=0.5, width=bar_width)
        #add text with the percentage of particles on the left of the bar
        ax[ind_m+1,2].text(0.8, 0.5, f"{round(psubducted*100, 1)}%", fontsize=10, ha="center", va="center", weight="bold", color = "mediumblue")
        ax[ind_m+1,2].text(0.8, psubducted + 0.5*pexhumed, f"{round(pexhumed*100, 1)}%", fontsize=10, ha="center", va="top", weight="bold", color = "indigo")
        ax[ind_m+1,2].text(0.8, psubducted + pexhumed + 0.5*pstagnant, f"{round(pstagnant*100, 1)}%", fontsize=10, ha="center", va="bottom", weight="bold", color = "darkgreen")
        

        ax[ind_m+1,2].axis("off")




    ax[0,2].axis("off")
    ax[n_models,0].set_ylabel("Duration (Myr)")
    ax[n_models,0].set_xlabel("Time (Myr)")
    ax[n_models,0].set_xlim(0, 50)
    ax[n_models,0].set_xticks([0, 10, 20, 30, 40, 50])

    ax[n_models,1].set_ylabel("Duration (Myr)")
    ax[n_models,1].set_xlabel("Time (Myr)")
    ax[n_models,1].set_xlim(0, 50)
    ax[n_models,1].set_xticks([0, 10, 20, 30, 40, 50])

    handles, labels = ax[1, 2].get_legend_handles_labels()

    # Place the legend in ax[0, 2]
    ax[0, 2].legend(handles, labels, loc='upper right', fontsize=10)


    # Add no vertical space between the plots
    plt.tight_layout()
    plt.subplots_adjust(hspace=0, wspace=0.1)

    plt.savefig(f"{figloc}/scatterplots_comparison.png", dpi=500)
    plt.close()





if __name__ == "__main__":
    main()
