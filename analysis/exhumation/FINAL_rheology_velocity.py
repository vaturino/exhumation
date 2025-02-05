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

    # colors
    color = [ "orangered", "firebrick"]


    cr_color = ["maroon", "orange"]



    fig, ax = plt.subplots(7, 2, figsize=(10, 7), height_ratios=[0.25, 0.1,  1, 1, 0.1, 1, 1])


    ax[1,0].set_visible(False)
    ax[1,1].set_visible(False)
    ax[4,0].set_visible(False)
    ax[4,1].set_visible(False)


    weak_idx = 0
    strong_idx = 0


    weak_stag = 5
    strong_stag = 5


    cr_cttV = pd.read_csv(f"/home/vturino/PhD/projects/exhumation/plots/single_models/kinematic_mu0.13_basalt7.5km_sed1km_SR1.8e-12_cttV/txt_files/2D_v.txt", sep="\s+")
    cr_Vto0 = pd.read_csv(f"/home/vturino/PhD/projects/exhumation/plots/single_models/kinematic_mu0.13_basalt7.5km_sed1km_SR1.8e-12_Vto0/txt_files/2D_v.txt", sep="\s+")

    cr_cttV.iloc[0] = np.nan
    cr_Vto0.iloc[0] = np.nan


    ax[0,0].plot(cr_Vto0["time"]/1e6, cr_Vto0["conv_rate"], color=cr_color[1], label=configs["names"][2])
    ax[0,0].plot(cr_cttV["time"]/1e6, cr_cttV["conv_rate"], color=cr_color[0], label=configs["names"][0])
    ax[0,0].axvline(x=35, color="grey", linestyle="--", linewidth = 1)
    ax[0,0].set_xlim(0, 50)
    ax[0,0].set_ylim(0, 8.5)
    ax[0,0].label_outer()

    ax[0,1].plot(cr_Vto0["time"]/1e6, cr_Vto0["conv_rate"], color=cr_color[1], label=configs["names"][2])
    ax[0,1].plot(cr_cttV["time"]/1e6, cr_cttV["conv_rate"], color=cr_color[0], label=configs["names"][0])
    ax[0,1].axvline(x=35, color="grey", linestyle="--", linewidth = 1)
    ax[0,1].set_xlim(0, 50)
    ax[0,1].set_ylim(0, 8.5)
    ax[0,1].label_outer()

    
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
                color=color[0],
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
                color=color[1],
                alpha=0.5,
                linewidth=1,
                element="step",
                legend=True,
            )
            ax[strong_idx+2, 1].axvline(x=35, color="grey", linestyle="--", linewidth = 1)
            ax[strong_idx+2, 1].set_yscale("log")
            ax[strong_idx+2, 1].set_xlim(0, 50)
            ax[strong_idx+2, 1].set_ylim(0.7, 5.e3)
            # ax[strong_idx+2, 1].set_ylabel("")
            # ax[strong_idx+2, 1].set_xticks([])
            # ax[strong_idx+2, 1].set_xlabel("")
            ax[strong_idx+2, 1].label_outer()
            # ax[strong_idx+2, 1].set_yticks([])
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

        # print(len(stagnant_list_expanded[stagnant_list_expanded["ti"]<24]))

        if configs["strength"][m] == "weak":
            sns.histplot(
                stagnant_list_expanded,
                x="ti",
                ax=ax[weak_stag, 0],
                bins=bin_edges,
                color=color[0],
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
                color=color[1],
                alpha=0.5,
                linewidth=1,
                element="step",
                legend=True
            )
            ax[strong_stag, 1].axvline(x=35, color="grey", linestyle="--", linewidth = 1)
            ax[strong_stag, 1].set_yscale("log")
            ax[strong_stag, 1].set_xlim(0, 50)
            ax[strong_stag, 1].set_ylim(0.7, 5.e3)
            # ax[strong_stag, 1].set_ylabel("")
            # ax[strong_stag, 1].set_xticks([])
            # ax[strong_stag, 1].set_xlabel("")
            # ax[strong_stag, 1].set_yticks([])
            ax[strong_stag, 1].label_outer()
            ax[strong_stag, 1].text(0.2, 0.8, configs["names"][m], horizontalalignment='center', verticalalignment='center', transform=ax[strong_stag, 1].transAxes)


            strong_stag += 1

    
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.0, wspace=0.0)
    plt.savefig(f"{plotloc}/vel_var_visc.pdf")
    plt.close()





if __name__ == "__main__":
    main()
