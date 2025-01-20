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
    linewidth = 1

 
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
    ax[0,0].set_yticks([0, 5])

    # # Histograms for tmax and tfin (Manual layering with plt.bar)
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


   

    # COLUMN 2

    # Plot 0: convergence rate
    ax[0,1].plot(cr["time"]/1e6, cr["conv_rate"], color="grey", label="Convergence rate (cm/yr)", linewidth=1)
    ax[0,1].axvline(x=35, color="grey", linestyle="--", linewidth = linewidth)
    ax[0,1].label_outer()  # Only show outer labels and tick labels
    ax[0,1].legend(loc='upper right')
    ax[0,1].set_ylim(0, cr["conv_rate"].max()+1)
    ax[0,1].set_xlim(0, 50)
    ax[0,1].set_title("Stagnant particles", fontsize=12)
    # ax[0,1].set_yticks([0, 5])
    ax[0,1].set_yticks([])
    ax[0,1].set_xticks([])
    ax[0,1].label_outer()  # Only show outer labels and tick labels



    # Histograms for tmax and tfin (Manual layering with plt.bar)
    stagnant_list_sorted = stagnant_list.sort_values(by='lithology', ascending=True)

    sns.histplot(stagnant_list_sorted, x="ti_kin", hue="lithology", ax=ax[1,1], palette=colors_tfin, alpha=alpha, linewidth=1, element="step", bins= bin_edges, legend = False)
    sns.histplot(stagnant_list_sorted, x="ti_dyn", hue="lithology", ax=ax[1,1], palette=colors_tfin, alpha=alpha, linewidth=1, element="step", bins= bin_edges, legend = False)
    sns.histplot(stagnant_list_sorted, x="ti_trans", hue="lithology", ax=ax[1,1], palette=colors_tfin, alpha=alpha, linewidth=1, element="step", bins= bin_edges, legend=True)

    ax[1,1].text(17, 2500, "Beginning of stagnation", fontsize=12)

    # Add a title, labels, and legend
    ax[1,1].set_xlabel("")
    ax[1,1].set_ylabel("")
    ax[1,1].set_xlim(0, 50)
    ax[1,1].set_ylim(0.7, 5.e3)
    ax[1,1].set_xticks([])
    ax[1,1].set_yticks([])
    ax[1,1].label_outer()  # Only show outer labels and tick labels
    #put y ticks to the right
    # ax[1,1].yaxis.tick_right()
    ax[1,1].set_yscale("log")
    ax[1,1].axvline(x=35, color="grey", linestyle="--", label="35 Ma", linewidth = linewidth)







    #print the particles in each bin of the histogram for exhumed particles
    new_bins = np.arange(0, 50, 5)
    offset = 0

    # for each lithology in the exhumed list, create a dataframe with the average values for each bin
    for l, lith in enumerate(exhumed_list["lithology"].unique()):
        lith_df = exhumed_list[exhumed_list["lithology"] == lith]
        avg_bins = pd.DataFrame(columns=["bin", "bin_edges", "particles_num", "avg_ti", "std_ti", "avg_duration", "std_duration"])
        current_bin = 1
        for i in range(0, len(new_bins)-1):
            avg_bins.loc[i, "bin_edges"] = f"{new_bins[i]}-{new_bins[i+1]}"
            avg_bins.loc[i, "particles_num"] = len(lith_df[(lith_df["ti"] >= new_bins[i]) & (lith_df["ti"] < new_bins[i+1])])
            avg_bins.loc[i, "avg_ti"] = lith_df[(lith_df["ti"] >= new_bins[i]) & (lith_df["ti"] < new_bins[i+1])]["ti"].mean()
            avg_bins.loc[i, "std_ti"] = lith_df[(lith_df["ti"] >= new_bins[i]) & (lith_df["ti"] < new_bins[i+1])]["ti"].std()
            avg_bins.loc[i, "avg_duration"] = lith_df[(lith_df["ti"] >= new_bins[i]) & (lith_df["ti"] < new_bins[i+1])]["time_interval"].mean()
            avg_bins.loc[i, "std_duration"] = lith_df[(lith_df["ti"] >= new_bins[i]) & (lith_df["ti"] < new_bins[i+1])]["time_interval"].std()
            if avg_bins.loc[i, "particles_num"] > 0:
                avg_bins.loc[i, "bin"] = current_bin + offset
                current_bin += 1
            else:
                avg_bins.loc[i, "bin"] = np.nan

        offset += len(avg_bins[avg_bins["particles_num"] > 0])
        ax[2,0].plot([avg_bins["avg_ti"], avg_bins["avg_ti"]+avg_bins["avg_duration"]], [avg_bins["bin"], avg_bins["bin"]], color=colors_tfin[lith], linewidth=1.5, linestyle='--', zorder = 1)
        ax[2,0].scatter(avg_bins["avg_ti"], avg_bins["bin"], color=colors_tfin[lith], label=lith, zorder = 10, marker=">", s = 100)
        ax[2,0].errorbar(
            avg_bins["avg_ti"] + avg_bins["avg_duration"],
            avg_bins["bin"], 
            xerr=avg_bins["std_duration"],
            fmt='o', 
            ecolor=colors_tin[lith], 
            color=colors_tin[lith],  
            markersize=5, 
            linewidth=1.5, 
            capsize=3
        )

    ax[2,0].invert_yaxis()
    ax[2,0].set_yticks([])
    ax[2,0].label_outer()  # Only show outer labels and tick labels
    ax[2,0].axvline(x=35, color="grey", linestyle="--", label="35 Ma", linewidth = linewidth)
    ax[2,0].set_xlim(0, 50)
    ax[2,0].set_xlabel("Time interval (Ma)")



    # do the same for stagnant particles: do not differentiate between ti_kin, ti_dyn, ti_trans etc
    new_bins = np.arange(0, 50, 5)
    offset = 0

    # for each lithology in the stagnant list, create a dataframe with the average values for each bin
    for l, lith in enumerate(stagnant_list["lithology"].unique()):
        lith_df = stagnant_list[stagnant_list["lithology"] == lith]
        avg_bins = pd.DataFrame(columns=["bin", "bin_edges", "particles_num", "avg_ti", "std_ti", "avg_duration", "std_duration"])
        current_bin = 1
        for i in range(0, len(new_bins)-1):
            avg_bins.loc[i, "bin_edges"] = f"{new_bins[i]}-{new_bins[i+1]}"
            particles_num = 0
            avg_ti = []
            avg_duration = []
            for ti_col, duration_col in [("ti_kin", "time_interval_kin"), ("ti_dyn", "time_interval_dyn"), ("ti_trans", "time_interval_trans")]:
                particles_num += len(lith_df[(lith_df[ti_col] >= new_bins[i]) & (lith_df[ti_col] < new_bins[i+1])])
                avg_ti.extend(lith_df[(lith_df[ti_col] >= new_bins[i]) & (lith_df[ti_col] < new_bins[i+1])][ti_col].tolist())
                avg_duration.extend(lith_df[(lith_df[ti_col] >= new_bins[i]) & (lith_df[ti_col] < new_bins[i+1])][duration_col].tolist())
            avg_bins.loc[i, "particles_num"] = particles_num
            avg_bins.loc[i, "avg_ti"] = np.mean(avg_ti) if avg_ti else np.nan
            avg_bins.loc[i, "std_ti"] = np.std(avg_ti) if avg_ti else np.nan
            avg_bins.loc[i, "avg_duration"] = np.mean(avg_duration) if avg_duration else np.nan
            avg_bins.loc[i, "std_duration"] = np.std(avg_duration) if avg_duration else np.nan
            if avg_bins.loc[i, "particles_num"] > 0:
                avg_bins.loc[i, "bin"] = current_bin + offset
                current_bin += 1
            else:
                avg_bins.loc[i, "bin"] = np.nan
            
        # print(avg_bins)

        offset += len(avg_bins[avg_bins["particles_num"] > 0])
        ax[2,1].plot([avg_bins["avg_ti"], avg_bins["avg_ti"]+avg_bins["avg_duration"]], [avg_bins["bin"], avg_bins["bin"]], color=colors_tfin[lith], linewidth=1.5, linestyle='--', zorder=2)
        ax[2,1].scatter(avg_bins["avg_ti"], avg_bins["bin"], color=colors_tfin[lith], label=lith, zorder=10, marker=">", s=100)
        ax[2,1].errorbar(
            avg_bins["avg_ti"] + avg_bins["avg_duration"],
            avg_bins["bin"], 
            xerr=avg_bins["std_duration"],
            fmt='o', 
            ecolor=colors_tin[lith], 
            color=colors_tin[lith],  
            markersize=5, 
            linewidth=1.5, 
            capsize=3,
            zorder=2
        )


    ax[2,1].invert_yaxis()
    ax[2,1].set_yticks([])
    ax[2,1].label_outer()  # Only show outer labels and tick labels
    ax[2,1].axvline(x=35, color="grey", linestyle="--", label="35 Ma", linewidth = linewidth)
    ax[2,1].set_xlim(0, 50)
    ax[2,1].set_xlabel("Time interval (Ma)")

   


    


    # Add no vertical space between the plots
    plt.tight_layout()
    plt.subplots_adjust(hspace=0, wspace=0.0)




    

    plt.savefig(f"{plot_loc}/timeline.png", dpi=500)
    plt.close()





if __name__ == "__main__":
    main()
