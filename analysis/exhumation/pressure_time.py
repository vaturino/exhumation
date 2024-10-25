#!/usr/bin/env python3

import pandas as pd
import numpy as np
import json
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import os

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
    stagnant_list = pd.read_csv(f"{txt_loc}/timing_stagnant_particles.txt", sep="\s+")
    cr = pd.read_csv(f"{txt_loc}/2D_v.txt", sep="\s+")
    cr["conv_rate"].iloc[0] = np.nan


    # Create a mapping for lithology colors
     # Define color palettes
    colors_tin = {
        "sed": "midnightblue",
        "oc": "saddlebrown",
        "ecl": "darkgreen",
        "serp": "maroon"
    }

    colors_tmax = {
        "sed": "mediumblue",
        "oc": "darkorange",
        "ecl": "forestgreen",
        "serp": "firebrick"
    }

    colors_tfin = {
        "sed": "cornflowerblue",
        "oc": "goldenrod",
        "ecl": "olivedrab",
        "serp": "indianred"
    }

    # Categorize 'vbur' and 'vexh' into 3 groups and capture bin edges
    # exhumed_list["vbur_category"], e_bin_edges_bur = pd.qcut(exhumed_list["vbur"], q=3, labels=["low", "mid", "high"], retbins=True)
    # e_bin_edges_bur = np.round(e_bin_edges_bur, 2)

    # exhumed_list["vexh_category"], e_bin_edges_exh = pd.qcut(exhumed_list['vexh'], q=3, labels=["low", "mid", "high"], retbins=True)
    # e_bin_edges_exh = np.round(e_bin_edges_exh, 2)

    # stagnant_list["vbur_category"], s_bin_edges_bur = pd.qcut(stagnant_list["vbur"], q=3, labels=["low", "mid", "high"], retbins=True)
    # s_bin_edges_bur = np.round(s_bin_edges_bur, 2)

    # stagnant_list["vstag_category"], s_bin_edges_stag = pd.qcut(stagnant_list['vstag'], q=3, labels=["low", "mid", "high"], retbins=True)
    # s_bin_edges_stag = np.round(s_bin_edges_stag, 2)


    # # Define marker sizes for each category
    # sizes = {"low": 20, "mid": 50, "high": 80}  # Adjust sizes for clarity
    # markers = {"low": "s", "mid": "X", "high": "d"}

    # Create the plot
    fig, ax = plt.subplots(2, 1, figsize=(15, 10))
    fig.suptitle("Pressure-time conditions")





    # Scatter plots for exhumed particles
    sns.scatterplot(data=exhumed_list, x="tin", y="Pin", hue="lithology", ax=ax[0], palette=colors_tin, legend=True)
    sns.scatterplot(data=exhumed_list, x="tmax", y="maxP", hue="lithology", ax=ax[0], palette=colors_tmax, legend=False)
    sns.scatterplot(data=exhumed_list, x="tfin", y="Pexh", hue="lithology", ax=ax[0], palette=colors_tfin, legend=False)
    # sns.scatterplot(data=exhumed_list, x="tin", y="maxP", hue="lithology", ax=ax[0], palette=colors_tin, legend=True, size="vbur_category", sizes=sizes)
    # sns.scatterplot(data=exhumed_list, x="tmax", y="maxP", hue="lithology", ax=ax[0], palette=colors_tmax, legend=False)
    # sns.scatterplot(data=exhumed_list, x="tfin", y="maxP", hue="lithology", ax=ax[0], palette=colors_tfin, legend=False, style="vexh_category", markers=markers)
    a1 = ax[0].twinx()
    a1.plot(cr["time"]/1e6, cr["conv_rate"], color="grey", label="Convergence rate", linewidth=3)
    a1.set_ylabel("Convergence rate (cm/yr)", color="grey", fontweight="bold")
    ax[0].set_xlabel("Time (Ma)")
    ax[0].set_ylabel("Pressure (GPa)")
    ax[0].set_title("Exhumed particles")
    ax[0].set_ylim(0, stagnant_list["maxP"].max()+0.2)

    # Create a custom legend for marker sizes
    handles, labels = ax[0].get_legend_handles_labels()
    unique_labels = set(labels)
    new_handles = []
    new_labels = []

    # Add lithology colors
    # for label in unique_labels:
    #     if label in colors_tin.keys():
    #         new_handles.append(plt.Line2D([0], [0], marker='o', color='w', label=label,
    #                                         markerfacecolor=colors_tin[label], markersize=10))
    #         new_labels.append(label)

    # # Add marker sizes for vbur categories with bin values
    # for category, size in sizes.items():
    #     low_value = e_bin_edges_bur[list(sizes.keys()).index(category)]
    #     high_value = e_bin_edges_bur[list(sizes.keys()).index(category) + 1]
    #     new_handles.append(plt.Line2D([0], [0], marker='o', color='w', label=f"vbur {category} ({low_value}, {high_value})",
    #                                     markerfacecolor='black', markersize=size / 7))  # Scale down for consistency
    #     new_labels.append(f"vbur {category} ({low_value}, {high_value})")

    # # Add marker sizes for vexh categories with bin values
    # for category, marker in markers.items():
    #     low_value = e_bin_edges_exh[list(markers.keys()).index(category)]
    #     high_value = e_bin_edges_exh[list(markers.keys()).index(category) + 1]
    #     new_handles.append(plt.Line2D([0], [0], marker=marker, color='w', label=f"vexh {category} ({low_value}, {high_value})",
    #                                     markerfacecolor='black', markersize=7))  # Scale down for consistency
    #     new_labels.append(f"vexh {category} ({low_value}, {high_value})")

    # ax[0].legend(handles=new_handles, labels=new_labels, title="Lithology and velocities", loc='upper right')



    # Scatter plots for stagnant particles
    # sns.scatterplot(data=stagnant_list, x="tin", y="maxP", hue="lithology", ax=ax[1], palette=colors_tin, legend=True, size="vbur_category", sizes=sizes)
    # sns.scatterplot(data=stagnant_list, x="tmax", y="maxP", hue="lithology", ax=ax[1], palette=colors_tmax, legend=False)
    # sns.scatterplot(data=stagnant_list, x="tfin", y="maxP", hue="lithology", ax=ax[1], palette=colors_tfin, legend=False, style="vstag_category", markers=markers)
    sns.scatterplot(data=stagnant_list, x="tin", y="Pin", hue="lithology", ax=ax[1], palette=colors_tin, legend=True)
    sns.scatterplot(data=stagnant_list, x="tmax", y="maxP", hue="lithology", ax=ax[1], palette=colors_tmax, legend=False)
    sns.scatterplot(data=stagnant_list, x="tfin", y="Pstag", hue="lithology", ax=ax[1], palette=colors_tfin, legend=False)
    a2 = ax[1].twinx()
    a2.plot(cr["time"]/1e6, cr["conv_rate"], color="grey", label="Convergence rate", linewidth=3)
    a2.set_ylabel("Convergence rate (cm/yr)", color="grey", fontweight="bold")
    ax[1].set_xlabel("Time (Ma)")
    ax[1].set_ylabel("Pressure (GPa)")
    ax[1].set_title("Stagnant particles")
    ax[1].set_ylim(0, stagnant_list["maxP"].max()+0.2)   
    ax[1].legend(loc='upper center')

    # Create a custom legend for marker sizes
    handles, labels = ax[1].get_legend_handles_labels()
    unique_labels = set(labels)
    new_handles = []
    new_labels = []

    # # Add lithology colors
    # for label in unique_labels:
    #     if label in colors_tin.keys():
    #         new_handles.append(plt.Line2D([0], [0], marker='o', color='w', label=label,
    #                                         markerfacecolor=colors_tin[label], markersize=10))
    #         new_labels.append(label)
    
    # # Add marker sizes for vbur categories with bin values
    # for category, size in sizes.items():
    #     low_value = s_bin_edges_bur[list(sizes.keys()).index(category)]
    #     high_value = s_bin_edges_bur[list(sizes.keys()).index(category) + 1]
    #     new_handles.append(plt.Line2D([0], [0], marker='o', color='w', label=f"vbur {category} ({low_value}, {high_value})",
    #                                     markerfacecolor='black', markersize=size / 7))  # Scale down for consistency
    #     new_labels.append(f"vbur {category} ({low_value}, {high_value})")
    
    # # Add marker sizes for vexh categories with bin values
    # for category, marker in markers.items():
    #     low_value = s_bin_edges_stag[list(markers.keys()).index(category)]
    #     high_value = s_bin_edges_stag[list(markers.keys()).index(category) + 1]
    #     new_handles.append(plt.Line2D([0], [0], marker=marker, color='w', label=f"vstag {category} ({low_value}, {high_value})",
    #                                     markerfacecolor='black', markersize=7))  # Scale down for consistency
    #     new_labels.append(f"vstag {category} ({low_value}, {high_value})")
    
    # ax[1].legend(handles=new_handles, labels=new_labels, title="Lithology and velocities", loc=(0.5, 0.37))

    

    plt.savefig(f"{plot_loc}/pressure_time.png", dpi=500)

if __name__ == "__main__":
    main()

