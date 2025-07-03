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

    xmax_plot = 65
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
    fig,ax = plt.subplots(6, 1, figsize=(10, 10), height_ratios=[0.25, 0.05, 1, 0.05, 1, 1])

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
    sns.histplot(exhumed_list_sorted, x="texh", hue="lithology", bins=bin_edges, ax=ax[2], palette=colors_tin, alpha=1, linewidth=1, element="step", fill=False, legend=False)
    

    # Set scale, labels, and limits
    ax[2].set_yscale("log")
    ax[2].label_outer()  # Only show outer labels and tick labels
    # ax[2].text(1, 20, " Particles \n Subduction", fontsize=18)
    ax[2].set_xlabel("")
    ax[2].set_ylabel("Particles count (log)")
    ax[2].set_xlim(0, xmax_plot)
    ax[2].set_ylim(0.7, 5.e3)
    ax[2].axvline(x=35, color="grey", linestyle="--", label="35 Ma", linewidth = linewidth)
    ax[2].text(0.32, -2, "Time at exhumable threshold", fontsize=12, transform=ax[1].transAxes)


    ax[3].set_visible(False)


    # # Histograms for tmax and tfin (Manual layering with plt.bar)
    # stagnant_list_sorted = stagnant_list.sort_values(by='lithology', ascending=True)

    # sns.histplot(stagnant_list_sorted, x="ti_kin", hue="lithology", ax=ax[4], palette=colors_tfin, alpha=alpha, linewidth=1, element="step", bins= bin_edges, legend = False)
    # sns.histplot(stagnant_list_sorted, x="ti_dyn", hue="lithology", ax=ax[4], palette=colors_tfin, alpha=alpha, linewidth=1, element="step", bins= bin_edges, legend = False)
    # sns.histplot(stagnant_list_sorted, x="ti_trans", hue="lithology", ax=ax[4], palette=colors_tfin, alpha=alpha, linewidth=1, element="step", bins= bin_edges, legend=True)

    # Initialize new columns
    stagnant_list["ti"] = np.nan
    stagnant_list["time_interval"] = np.nan
    stagnant_list["lith_time"] = np.nan

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

    # print(stagnant_list[(stagnant_list["lithology"] == "oc") & (stagnant_list["ti"] <18)])


    sns.histplot(stagnant_list_expanded, x="ti", ax=ax[4], hue = "lith_time", palette= colors_tfin, alpha=alpha, linewidth=1, element="step", bins= bin_edges, legend=False)
    
    
    ax[4].text(17, 2400, "Beginning of stagnation", fontsize=12)

    # Add a title, labels, and legend
    ax[4].set_xlabel("")
    ax[4].set_ylabel("")
    ax[4].set_xlim(0, xmax_plot)
    ax[4].set_ylim(0.7, 5.e3)
    ax[4].set_xticks([])
    ax[4].set_yticks([])
    ax[4].label_outer()  # Only show outer labels and tick labels
    #put y ticks to the right
    # ax[4].yaxis.tick_right()
    ax[4].set_yscale("log")
    ax[4].axvline(x=35, color="grey", linestyle="--", label="35 Ma", linewidth = linewidth)







    new_bins = np.arange(0, time_max, 5)
    offset = 0


    # for each lithology in the stagnant list, create a dataframe with the average values for each bin
    # for l, lith in enumerate(stagnant_list["lithology"].unique()):
    for l, lith in enumerate(["sed", "oc", "ecl"]):
        # print(f"Processing lithology: {lith}")
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
        ax[5].plot([avg_bins["avg_ti"], avg_bins["avg_ti"]+avg_bins["avg_duration"]], [avg_bins["bin"], avg_bins["bin"]], color=colors_tfin[lith], linewidth=1.5, linestyle='--', zorder=2)
        ax[5].scatter(avg_bins["avg_ti"], avg_bins["bin"], color=colors_tfin[lith], label=lith, zorder=10, marker=">", s=100, linewidth = 0)
        ax[5].errorbar(
            avg_bins["avg_ti"] + avg_bins["avg_duration"],
            avg_bins["bin"], 
            xerr=avg_bins["std_duration"],
            fmt='o', 
            ecolor=colors_tfin[lith], 
            color=colors_tfin[lith],  
            markersize=7, 
            linewidth=0, 
            # capsize=5,
            # capthick=2,
            zorder=2
        )


    ax[5].invert_yaxis()
    ax[5].set_yticks([])
    ax[5].label_outer()  # Only show outer labels and tick labels
    ax[5].axvline(x=35, color="grey", linestyle="--", label="35 Ma", linewidth = linewidth)
    ax[5].set_xlim(0, xmax_plot)
    ax[5].set_xlabel("Time (Myr)")

   
    plt.subplots_adjust(hspace=0.0)





    

    plt.savefig(f"{plot_loc}/timeline.pdf", dpi=500)
    plt.close()





if __name__ == "__main__":
    main()
