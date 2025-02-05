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
    color = ["mediumslateblue", "darkslateblue", "orangered", "firebrick", "forestgreen", "darkgreen"]

    colors_tfin = {
        "sed": "cornflowerblue",
        "oc": "#E3B64F",
        "ecl": "#A0C93D",
        "serp": "lightsalmon"
    }


    # fig, ax = plt.subplots(9, 2, figsize=(10, 7), height_ratios=[0.25, 0.1,  1, 1, 1, 0.1, 1, 1, 1])
    
    # # Plot 0: convergence rate
    # ax[0,0].plot(cr["time"]/1e6, cr["conv_rate"], color="grey", label="Convergence rate (cm/yr)", linewidth=1)
    # ax[0,0].axvline(x=35, color="grey", linestyle="--", linewidth = 1)
    # ax[0,0].label_outer()  # Only show outer labels and tick labels
    # ax[0,0].set_xlim(0, 50)

    # ax[0,1].plot(cr["time"]/1e6, cr["conv_rate"], color="grey", label="Convergence rate (cm/yr)", linewidth=1)
    # ax[0,1].axvline(x=35, color="grey", linestyle="--", linewidth = 1)
    # ax[0,1].label_outer()  # Only show outer labels and tick labels
    # ax[0,1].set_xlim(0, 50)

    # ax[1,0].set_visible(False)
    # ax[1,1].set_visible(False)
    # ax[5,0].set_visible(False)
    # ax[5,1].set_visible(False)





    # fig.subplots_adjust(hspace=0.0, wspace=0.0)

    # weak_idx = 0
    # strong_idx = 0


    # weak_stag = 6
    # strong_stag = 6


    # for m, mod in enumerate(configs['models']):
    #     mod_loc = f"/home/vturino/PhD/projects/exhumation/plots/single_models/{mod}"
    #     exhumed_list = pd.read_csv(f"{mod_loc}/txt_files/exhumed_times.txt", sep="\s+")
    #     stagnant_list = pd.read_csv(f"{mod_loc}/txt_files/stagnant_particles.txt", sep="\s+")

    #     bin_edges = np.arange(0, 50, 1)  # Define bin edges

    #     if configs["strength"][m] == "weak":
    #         sns.histplot(
    #             exhumed_list,
    #             x="ti",
    #             ax=ax[weak_idx+2, 0],
    #             bins=bin_edges,
    #             hue = "lithology",
    #             palette = colors_tfin,
    #             alpha=0.5,
    #             linewidth=1,
    #             element="step",
    #             legend=True
    #         )

    #         ax[weak_idx+2, 0].axvline(x=35, color="grey", linestyle="--", linewidth = 1)
    #         ax[weak_idx+2, 0].set_yscale("log")
    #         ax[weak_idx+2, 0].set_xlim(0, 50)
    #         ax[weak_idx+2, 0].set_ylim(0.7, 5.e3)
    #         ax[weak_idx+2, 0].set_ylabel("Count (log)")
    #         ax[weak_idx+2, 0].set_xticks([])
    #         ax[weak_idx+2, 0].set_xlabel("")
    #         ax[weak_idx+2, 0].text(0.2, 0.8, configs["names"][m], horizontalalignment='center', verticalalignment='center', transform=ax[weak_idx+2, 0].transAxes)


    #         weak_idx += 1

    #     else:

    #         sns.histplot(
    #             exhumed_list,
    #             x="ti",
    #             ax=ax[strong_idx+2, 1],
    #             bins=bin_edges,
    #             hue = "lithology",
    #             palette = colors_tfin,
    #             alpha=0.5,
    #             linewidth=1,
    #             element="step",
    #             legend=True,
    #         )
    #         ax[strong_idx+2, 1].axvline(x=35, color="grey", linestyle="--", linewidth = 1)
    #         ax[strong_idx+2, 1].set_yscale("log")
    #         ax[strong_idx+2, 1].set_xlim(0, 50)
    #         ax[strong_idx+2, 1].set_ylim(0.7, 5.e3)
    #         ax[strong_idx+2, 1].set_ylabel("")
    #         ax[strong_idx+2, 1].set_xticks([])
    #         ax[strong_idx+2, 1].set_xlabel("")
    #         ax[strong_idx+2, 1].set_yticks([])
    #         ax[strong_idx+2, 1].text(0.2, 0.8, configs["names"][m], horizontalalignment='center', verticalalignment='center', transform=ax[strong_idx+2, 1].transAxes)

    #         strong_idx += 1

    #     #stagnant particles
    #     # Handle ti_kin, ti_dyn, ti_trans
    #     new_rows = []
    #     for index, row in stagnant_list.iterrows():
    #         ti_values = [row["ti_dyn"], row["ti_kin"], row["ti_trans"]]
    #         time_interval_values = [row["time_interval_dyn"], row["time_interval_kin"], row["time_interval_trans"]]
    #         lithology_values = [row["lithology_dyn"], row["lithology_kin"], row["lithology_trans"]]
            
    #         for ti, time_interval, lithology in zip(ti_values, time_interval_values, lithology_values):
    #             if not np.isnan(ti):  # Only include non-NaN values
    #                 new_row = row.copy()
    #                 new_row["ti"] = ti
    #                 new_row["time_interval"] = time_interval
    #                 new_row["lith_time"] = lithology
    #                 new_rows.append(new_row)

    #     # Concatenate new rows to create the final DataFrame
    #     stagnant_list_expanded = pd.DataFrame(new_rows)

    #     if configs["strength"][m] == "weak":
    #         sns.histplot(
    #             stagnant_list_expanded,
    #             x="ti",
    #             ax=ax[weak_stag, 0],
    #             bins=bin_edges,
    #             hue="lith_time",
    #             palette=colors_tfin,
    #             alpha=0.5,
    #             linewidth=1,
    #             element="step",
    #             legend=True,
    #         )

    #         ax[weak_stag, 0].axvline(x=35, color="grey", linestyle="--", linewidth = 1)
    #         ax[weak_stag, 0].set_yscale("log")
    #         ax[weak_stag, 0].set_xlim(0, 50)
    #         ax[weak_stag, 0].set_ylim(0.7, 5.e3)
    #         ax[weak_stag, 0].set_ylabel("Count (log)")
    #         ax[weak_stag, 0].set_xticks([])
    #         ax[weak_stag, 0].set_xlabel("")
    #         ax[weak_stag, 0].text(0.2, 0.8, configs["names"][m], horizontalalignment='center', verticalalignment='center', transform=ax[weak_stag, 0].transAxes)


    #         weak_stag += 1

    #     else:
            
    #         sns.histplot(
    #             stagnant_list_expanded,
    #             x="ti",
    #             ax=ax[strong_stag, 1],
    #             bins=bin_edges,
    #             hue="lith_time",
    #             palette=colors_tfin,
    #             alpha=0.5,
    #             linewidth=1,
    #             element="step",
    #             legend=True
    #         )
    #         ax[strong_stag, 1].axvline(x=35, color="grey", linestyle="--", linewidth = 1)
    #         ax[strong_stag, 1].set_yscale("log")
    #         ax[strong_stag, 1].set_xlim(0, 50)
    #         ax[strong_stag, 1].set_ylim(0.7, 5.e3)
    #         ax[strong_stag, 1].set_ylabel("")
    #         ax[strong_stag, 1].set_xticks([])
    #         ax[strong_stag, 1].set_xlabel("")
    #         ax[strong_stag, 1].set_yticks([])
    #         ax[strong_stag, 1].text(0.2, 0.8, configs["names"][m], horizontalalignment='center', verticalalignment='center', transform=ax[strong_stag, 1].transAxes)


    #         strong_stag += 1




    #     ax[4,1].set_visible(False)
    #     ax[8,1].set_visible(False)
    #     # ax[8,0].legend(loc="best")



    
    # plt.tight_layout()
    # plt.subplots_adjust(hspace=0.0, wspace=0.0)
    # plt.savefig(f"{plotloc}/{plotname}")
    # plt.close()








    fig, ax = plt.subplots(9, 2, figsize=(20, 15), height_ratios=[0.3, 0.3,  1, 1, 1, 0.3, 1, 1, 1])
    
    # Plot 0: convergence rate
    ax[0,0].plot(cr["time"]/1e6, cr["conv_rate"], color="grey", label="Convergence rate (cm/yr)", linewidth=1)
    ax[0,0].axvline(x=35, color="grey", linestyle="--", linewidth = 1)
    ax[0,0].label_outer()  # Only show outer labels and tick labels
    ax[0,0].set_xlim(0, 50)
    ax[0,0].set_ylim(0, cr["conv_rate"].max()+1)

    ax[0,1].plot(cr["time"]/1e6, cr["conv_rate"], color="grey", label="Convergence rate (cm/yr)", linewidth=1)
    ax[0,1].axvline(x=35, color="grey", linestyle="--", linewidth = 1)
    ax[0,1].label_outer()  # Only show outer labels and tick labels
    ax[0,1].set_xlim(0, 50)
    ax[0,1].set_ylim(0, cr["conv_rate"].max()+1)

    ax[1,0].set_visible(False)
    ax[1,1].set_visible(False)
    ax[5,0].set_visible(False)
    ax[5,1].set_visible(False)






    weak_idx = 0
    strong_idx = 0


    weak_stag = 6
    strong_stag = 6


    for m, mod in enumerate(configs['models']):
        mod_loc = f"/home/vturino/PhD/projects/exhumation/plots/single_models/{mod}"
        exhumed_list = pd.read_csv(f"{mod_loc}/txt_files/exhumed_times.txt", sep="\s+")
        stagnant_list = pd.read_csv(f"{mod_loc}/txt_files/stagnant_particles.txt", sep="\s+")

        bin_edges = np.arange(0, 50, 1)  # Define bin edges

        if configs["strength"][m] == "weak":
            sns.histplot(
                exhumed_list,
                x="ti",
                ax=ax[weak_idx+2, 0],
                bins=bin_edges,
                color=color[m],
                alpha=0.5,
                linewidth=1,
                element="step",
                legend=True
            )

            ax[weak_idx+2, 0].axvline(x=35, color="grey", linestyle="--", linewidth = 1)
            ax[weak_idx+2, 0].set_yscale("log")
            ax[weak_idx+2, 0].set_xlim(0, 50)
            ax[weak_idx+2, 0].set_ylim(0.7, 5.e3)
            ax[weak_idx+2, 0].set_ylabel("Count (log)")
            ax[weak_idx+2, 0].set_xticks([])
            ax[weak_idx+2, 0].set_xlabel("")
            ax[weak_idx+2, 0].text(0.2, 0.8, configs["names"][m], horizontalalignment='center', verticalalignment='center', transform=ax[weak_idx+2, 0].transAxes)


            weak_idx += 1

        else:

            sns.histplot(
                exhumed_list,
                x="ti",
                ax=ax[strong_idx+2, 1],
                bins=bin_edges,
                color=color[m],
                alpha=0.5,
                linewidth=1,
                element="step",
                legend=True,
            )
            ax[strong_idx+2, 1].axvline(x=35, color="grey", linestyle="--", linewidth = 1)
            ax[strong_idx+2, 1].set_yscale("log")
            ax[strong_idx+2, 1].set_xlim(0, 50)
            ax[strong_idx+2, 1].set_ylim(0.7, 5.e3)
            ax[strong_idx+2, 1].set_ylabel("")
            ax[strong_idx+2, 1].label_outer()
            ax[strong_idx+2, 1].text(0.2, 0.8, configs["names"][m], horizontalalignment='center', verticalalignment='center', transform=ax[strong_idx+2, 1].transAxes)

            strong_idx += 1

        #stagnant particles
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

        if configs["strength"][m] == "weak":
            sns.histplot(
                stagnant_list_expanded,
                x="ti",
                ax=ax[weak_stag, 0],
                bins=bin_edges,
                color=color[m],
                alpha=0.5,
                linewidth=1,
                element="step",
                legend=True,
            )

            ax[weak_stag, 0].axvline(x=35, color="grey", linestyle="--", linewidth = 1)
            ax[weak_stag, 0].set_yscale("log")
            ax[weak_stag, 0].set_xlim(0, 50)
            ax[weak_stag, 0].set_ylim(0.7, 5.e3)
            ax[weak_stag, 0].set_ylabel("Count (log)")
            ax[weak_stag, 0].set_xticks([])
            ax[weak_stag, 0].set_xlabel("")
            ax[weak_stag, 0].text(0.2, 0.8, configs["names"][m], horizontalalignment='center', verticalalignment='center', transform=ax[weak_stag, 0].transAxes)


            weak_stag += 1

        else:
            
            sns.histplot(
                stagnant_list_expanded,
                x="ti",
                ax=ax[strong_stag, 1],
                bins=bin_edges,
                color=color[m],
                alpha=0.5,
                linewidth=1,
                element="step",
                legend=True
            )
            ax[strong_stag, 1].axvline(x=35, color="grey", linestyle="--", linewidth = 1)
            ax[strong_stag, 1].set_yscale("log")
            ax[strong_stag, 1].set_xlim(0, 50)
            ax[strong_stag, 1].set_ylim(0.7, 5.e3)
            ax[strong_stag, 1].set_ylabel("")
            ax[strong_stag, 1].label_outer()
            ax[strong_stag, 1].text(0.2, 0.8, configs["names"][m], horizontalalignment='center', verticalalignment='center', transform=ax[strong_stag, 1].transAxes)


            strong_stag += 1




        ax[4,1].set_visible(False)
        ax[8,1].set_visible(False)
        # ax[8,0].legend(loc="best")

        ax[8,0].set_xlabel("Time (Myr)")
        ax[8,0].set_xticks([0, 10, 20, 30, 40, 50])
        ax[7,1].set_xlabel("Time (Myr)")
        ax[7,1].set_xticks([0, 10, 20, 30, 40, 50])
        ax[3,1].set_xlabel("Time (Myr)")
        ax[3,1].set_xticks([0, 10, 20, 30, 40, 50])


    
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.0, wspace=0.1)
    plt.savefig(f"{plotloc}/rheo_no_lith.pdf")
    plt.close()





if __name__ == "__main__":
    main()
