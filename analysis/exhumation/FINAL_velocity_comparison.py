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
        # sns.histplot(exhumed_list_sorted, x="tm", hue="lithology", bins=bin_edges, ax=ax[ind_m+1,0], palette=colors_tfin, alpha=alpha, linewidth=1, element="step", legend = False)
        

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



        # # Histograms for tmax and tfin (Manual layering with plt.bar)
        # stagnant_list["ti"] = np.nan
        # stagnant_list["time_interval"] = np.nan
        # stagnant_list["lith_time"] = np.nan
        # # sns.histplot(stagnant_list, x="ti_dyn", ax=ax[ind_m+1,1], hue = "lithology", palette= colors_tfin, alpha=alpha, linewidth=1, element="step", bins= bin_edges, legend=False)

        # # Handle ti_kin, ti_dyn, ti_trans
        # new_rows = []
        # for index, row in stagnant_list.iterrows():
        #     ti_values = [row["ti_kin"], row["ti_dyn"], row["ti_trans"]]
        #     non_nan_ti_values = [ti for ti in ti_values if not np.isnan(ti)]
        #     if len(non_nan_ti_values) == 1:
        #         stagnant_list.at[index, "ti"] = non_nan_ti_values[0]
        #     elif len(non_nan_ti_values) > 1:
        #         for ti in non_nan_ti_values:
        #             new_row = row.copy()
        #             new_row["ti"] = ti
        #             new_rows.append(new_row)
        #         stagnant_list.at[index, "ti"] = np.nan  # Mark original row as NaN to be dropped later

        # # Concatenate new rows to the DataFrame and drop original rows with NaN "ti"
        # stagnant_list = pd.concat([stagnant_list, pd.DataFrame(new_rows)], ignore_index=True)
        # stagnant_list = stagnant_list.dropna(subset=["ti"])

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

        print(stagnant_list_expanded[(stagnant_list_expanded["lithology"] == "oc") & (stagnant_list_expanded["ti"] >22) & (stagnant_list_expanded["ti"] < 23)]["id"])


        sns.histplot(stagnant_list_expanded, x="ti", ax=ax[ind_m+1,1], hue = "lith_time", palette= colors_tfin, alpha=alpha, linewidth=1, element="step", bins= bin_edges, legend=False)

        ax[ind_m+1,1].text(17, 2500, "Beginning of stagnation", fontsize=12)

        # Add a title, labels, and legend
        ax[ind_m+1,1].set_xlabel("")
        ax[ind_m+1,1].set_ylabel("")
        ax[ind_m+1,1].set_xlim(13, 50)
        ax[ind_m+1,1].set_ylim(0.7, 5.e3)
        ax[ind_m+1,1].set_xticks([])
        ax[ind_m+1,1].set_yticks([])
        ax[ind_m+1,1].label_outer()  # Only show outer labels and tick labels
        #put y ticks to the right
        # ax[ind_m+1,1].yaxis.tick_right()
        ax[ind_m+1,1].set_yscale("log")
        ax[ind_m+1,1].axvline(x=35, color="grey", linestyle="--", label="35 Ma", linewidth = linewidth)


        # add bars with percentage of exhumed/stagnant particles

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


        bar_width = 0.02
        ax[ind_m+1,2].bar(0.8, psubducted, color=classification["subducted"], 
                        label="Subducted", edgecolor="black", linewidth=0.5, width=bar_width)
        ax[ind_m+1,2].bar(0.8, pexhumed, bottom=psubducted, color=classification["exhumed"], 
                        label="Exhumed", edgecolor="black", linewidth=0.5, width=bar_width)
        ax[ind_m+1,2].bar(0.8, pstagnant, bottom=psubducted+pexhumed, color=classification["stagnant"], 
                        label="Stagnant", edgecolor="black", linewidth=0.5, width=bar_width)
        #add text with the percentage of particles on the left of the bar
        ax[ind_m+1,2].text(0.8, 0.5, f"{round(psubducted*100, 0)}%", fontsize=10, ha="center", va="center", weight="bold", color = "mediumblue")
        ax[ind_m+1,2].text(0.8, psubducted + 0.5*pexhumed, f"{round(pexhumed*100, 0)}%", fontsize=10, ha="center", va="top", weight="bold", color = "indigo")
        ax[ind_m+1,2].text(0.8, psubducted + pexhumed + 0.5*pstagnant, f"{round(pstagnant*100, 0)}%", fontsize=10, ha="center", va="bottom", weight="bold", color = "darkgreen")
        

        ax[ind_m+1,2].axis("off")

    ax[0,2].axis("off")
    ax[count,0].set_xticks([20, 30, 40, 50])
    ax[count,0].set_xlabel("Time (Myr)")

    ax[count,1].set_xticks([20, 30, 40, 50])
    ax[count,1].set_xlabel("Time (Myr)")


    # Add no vertical space between the plots
    plt.tight_layout()
    plt.subplots_adjust(hspace=0, wspace=0.0)




    

    plt.savefig(f"{figloc}/{plotname}", dpi=500)
    plt.close()





if __name__ == "__main__":
    main()
