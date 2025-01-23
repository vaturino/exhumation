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

    cr_colors = [
        "maroon",
        "red",
        "orange"
    ]

    classification = {
        "subducted": "lightcyan",
        "exhumed": "orchid",
        "stagnant": "olivedrab"
    }


    alpha = 0.5
    linewidth = 1

    figloc = f"/home/vturino/PhD/projects/exhumation/plots/comparison/"
    plotname = f"velocity_comparison_timing_{configs['test']}.pdf"
    count = len(configs["models"])

    perc = pd.DataFrame(columns=["model", "lithology", "num", "perc", "rel_perc"])


    for ind_m, m in enumerate(configs["models"]):
        # Create the folders to save the plots
        plot_loc = f"/home/vturino/PhD/projects/exhumation/plots/single_models/{m}"
        txt_loc = f'{plot_loc}/txt_files'
        if not os.path.exists(txt_loc):
            os.mkdir(txt_loc)

        # Load the data
        exhumed_list = pd.read_csv(f"{txt_loc}/exhumed_times.txt", sep="\s+")
        stagnant_list = pd.read_csv(f"{txt_loc}/stagnant_particles.txt", sep="\s+")
        particles = pd.read_csv(f"{txt_loc}/particles_indexes.csv", sep="\s+")
        nparts = len(particles)

        # count amount of stagnant particles by lithology
        for lith in stagnant_list["lithology"].unique():
            num = len(stagnant_list[stagnant_list["lithology"] == lith])
            perc.loc[len(perc)] = [m, lith, num, num/nparts*100, 0]
    
    # Calculate relative percentage for each lithology
    for lith in perc["lithology"].unique():
        if ind_m == 0:
            perc.loc[perc["lithology"] == lith, "rel_perc"] = 100
        else:
            base_num = perc[(perc["model"] == configs["models"][0]) & (perc["lithology"] == lith)]["num"].values[0]
            perc.loc[perc["lithology"] == lith, "rel_perc"] = (perc[perc["lithology"] == lith]["num"] / base_num) * 100

    # Exaggerate relative percentage for each lithology
    for lith in perc["lithology"].unique():
        for ind_m, m in enumerate(configs["models"]):
            if ind_m > 0:
                if perc[(perc["model"] == m) & (perc["lithology"] == lith)]["rel_perc"].values[0] < 100:
                    perc.loc[(perc["model"] == m) & (perc["lithology"] == lith), "rel_perc"] = perc[(perc["lithology"] == lith)]["rel_perc"] /(1.2*ind_m)
                else:
                    perc.loc[(perc["model"] == m) & (perc["lithology"] == lith), "rel_perc"] = perc[(perc["model"] == m) & (perc["lithology"] == lith)]["rel_perc"] * 1.2*ind_m


    # # Create the plot
    fig,ax = plt.subplots(count+1, 3, figsize=(10, 6), height_ratios=[0.25]+[1]*count, width_ratios=[1,1,0.1])


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
        ax[0,0].plot(cr["time"]/1e6, cr["conv_rate"], color=cr_colors[ind_m], linewidth=1)
        ax[0,0].axvline(x=35, color='grey', linestyle="--", linewidth = linewidth)
        ax[0,0].label_outer()  # Only show outer labels and tick labels
        # ax[0,0].legend(loc='upper right')
        ax[0,0].set_ylim(0, cr["conv_rate"].max()+1)
        ax[0,0].set_xlim(13, 50)
        ax[0,0].set_title("Exhumable particles", fontsize=12)
        ax[0,0].set_yticks([0, 5])
        ax[0,0].set_xticks([])

        # # Histograms for tmax and tfin (Manual layering with plt.bar)
        bin_edges = np.arange(0, 50,1)  # Define bin edges


        # Sort exhumed_list by lithology to ensure 'oc' is plotted over 'sed'
        exhumed_list_sorted = exhumed_list.sort_values(by='lithology', ascending=True)
        sns.histplot(exhumed_list_sorted, x="ti", hue="lithology", bins=bin_edges, ax=ax[ind_m+1,0], palette=colors_tfin, alpha=alpha, linewidth=1, element="step", legend=False)
        

        # Set scale, labels, and limits
        ax[ind_m+1,0].set_yscale("log")
        ax[ind_m+1,0].label_outer()  # Only show outer labels and tick labels
        # ax[ind_m+1,0].text(1, 20, " Particles \n Subduction", fontsize=18)
        ax[ind_m+1,0].set_xlabel("")
        ax[ind_m+1,0].set_ylabel("Particles count (log)")
        ax[ind_m+1,0].set_xlim(13, 50)
        ax[ind_m+1,0].set_xticks([])
        ax[ind_m+1,0].set_ylim(0.7, 5.e3)
        ax[ind_m+1,0].axvline(x=35, color="grey", linestyle="--", label="35 Ma", linewidth = linewidth)
        ax[ind_m+1,0].text(0.32, 0.9, "Time at exhumable threshold", fontsize=12, transform=ax[ind_m+1,0].transAxes)


    

        # COLUMN 2

        # Plot 0: convergence rate
        ax[0,1].plot(cr["time"]/1e6, cr["conv_rate"], color= cr_colors[ind_m], linewidth=1)
        ax[0,1].axvline(x=35, color="grey", linestyle="--", linewidth = linewidth)
        ax[0,1].label_outer()  # Only show outer labels and tick labels
        # ax[0,1].legend(loc='upper right')
        ax[0,1].set_ylim(0, cr["conv_rate"].max()+1)
        ax[0,1].set_xlim(13, 50)
        ax[0,1].set_title("Stagnant particles", fontsize=12)
        # ax[0,1].set_yticks([0, 5])
        ax[0,1].set_yticks([])
        ax[0,1].set_xticks([])
        ax[0,1].label_outer()  # Only show outer labels and tick labels



        # duration of stagnation
        new_bins = np.arange(0, 50, 10)
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

            size = perc[(perc["model"] == m) & (perc["lithology"] == lith)]["rel_perc"].values[0]


            offset += len(avg_bins[avg_bins["particles_num"] > 0])
            ax[ind_m+1,1].plot([avg_bins["avg_ti"], avg_bins["avg_ti"]+avg_bins["avg_duration"]], [avg_bins["bin"], avg_bins["bin"]], color=colors_tfin[lith], linewidth=1.5, linestyle='--', zorder=2)
            ax[ind_m+1,1].scatter(avg_bins["avg_ti"], avg_bins["bin"], color=colors_tfin[lith], label=lith, zorder=10, marker=">", linewidth = 0, s=size)
            # sns.scatterplot(ax=ax[ind_m+1,1], data=avg_bins, x="avg_ti", y="bin", color=colors_tfin[lith], zorder=10, size="particles_num", legend=False, style="bin", markers=">", linewidth=0)
            ax[ind_m+1,1].errorbar(
                avg_bins["avg_ti"] + avg_bins["avg_duration"],
                avg_bins["bin"], 
                xerr=avg_bins["std_duration"],
                fmt='o', 
                ecolor=colors_tfin[lith], 
                color=colors_tfin[lith],  
                markersize=7, 
                linewidth=0, 
                capsize=5,
                capthick=2,
                zorder=2
            )


        ax[ind_m+1,1].invert_yaxis()
        ax[ind_m+1,1].set_yticks([])
        ax[ind_m+1,1].label_outer()  # Only show outer labels and tick labels
        ax[ind_m+1,1].axvline(x=35, color="grey", linestyle="--", linewidth = linewidth)
        ax[ind_m+1,1].set_xlim(13, 50)
        ax[ind_m+1,1].set_xlabel("Time (Myr)")








        # add bars with percentage of exhumed/stagnant particles




        bar_width = 0.02
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
    ax[count,0].set_xticks([20, 30, 40, 50])
    ax[count,0].set_xlabel("Time (Myr)")

    ax[count,1].set_xticks([20, 30, 40, 50])
    ax[count,1].set_xlabel("Time (Myr)")
    # ax[count,1].legend(loc='center left', fontsize=10)


    # Add no vertical space between the plots
    plt.tight_layout()
    plt.subplots_adjust(hspace=0, wspace=0.0)




    

    plt.savefig(f"{figloc}/{plotname}", dpi=500)
    plt.close()





if __name__ == "__main__":
    main()
