#!/usr/bin/env python3

import pandas as pd
import numpy as np
import json
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import os
import matplotlib.patches as mpatches


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


    labels_font = 13

    # Create the plot
    # fig, ax = plt.subplots(2, 1, figsize=(15, 10))
    #mosaic: the plot on top is thinner than the one on the bottom
    fig, ax = plt.subplots(2, 1, figsize=(15, 10), gridspec_kw={'height_ratios': [1, 2]})
    fig.suptitle("Pressure-time conditions")

    # Plot 0: convergence rate
    ax[0].plot(cr["time"]/1e6, cr["conv_rate"], color="grey", label="Convergence rate", linewidth=3)
    ax[0].set_ylabel("Convergence rate (cm/yr)", color="black", fontsize=labels_font)
    ax[0].set_xlabel("Time (Ma)", fontsize=labels_font)
    ax[0].set_title("Convergence rate and exhumed particles", fontsize=15)
    ax[0].set_ylim(0, cr["conv_rate"].max()+0.2)
    ax[0].set_xlim(0, 51)
    ax[0].label_outer()  # Only show outer labels and tick labels
    # vertical line at 35 Ma
    ax[0].axvline(x=35, color="grey", linestyle="--", label="35 Ma")

    # Plot 1: Pressure-time for exhumed particles
    ax[1].add_patch(plt.Rectangle((0, 0), exhumed_list["tin"].max(), 2000, color='grey', alpha=0.5))
    # Create legend handles manually for tfin and tmax using unique lithologies and colors
    unique_lithologies = exhumed_list["lithology"].unique()
    handles_tmax = [mpatches.Patch(color=colors_tmax[lith], label=lith) for lith in unique_lithologies]
    handles_tfin = [mpatches.Patch(color=colors_tfin[lith], label=lith) for lith in unique_lithologies]

    sns.histplot(
        data=exhumed_list, x="tmax", hue="lithology", ax=ax[1],
        palette=colors_tmax, legend=False, alpha=1, multiple="stack",
        edgecolor="black", linewidth=0.5
    )
    sns.histplot(
        data=exhumed_list, x="tfin", hue="lithology", ax=ax[1],
        palette=colors_tfin, legend=False, alpha=1, multiple="stack",
        edgecolor="black", linewidth=0.5
    )
    # Add custom legends for tmax and tfin
    legend_tmax = ax[1].legend(handles=handles_tmax, title="Peak Pressure", loc=(0.75, 0.85), frameon=True)
    ax[1].add_artist(legend_tmax)  # Manually add the tmax legend

    ax[1].legend(handles=handles_tfin, title="Time at Exhumability", loc=(0.85, 0.85), frameon=True)

    # Set scale, labels, and limits
    ax[1].set_yscale("log")
    ax[1].text(1, 20, " Particles \n Subduction", fontsize=18)
    ax[1].set_xlabel("Time (Ma)", fontsize=labels_font)
    ax[1].set_ylabel("Particles count (log)", fontsize=labels_font)
    ax[1].set_xlim(0, 51)
    ax[1].set_ylim(1, 2.e3)
    ax[1].axvline(x=35, color="grey", linestyle="--", label="35 Ma")

    # Add no vertical space between the plots
    plt.subplots_adjust(hspace=0)


    plt.savefig(f"{plot_loc}/exhumVScr_hist.png", dpi=500)
    plt.close()


    ################ TOTAL PLOT ################
    labels_font = 15

    # Create the plot
    # fig, ax = plt.subplots(2, 1, figsize=(15, 10))
    #mosaic: the plot on top is thinner than the one on the bottom
    fig, ax = plt.subplots(3, 1, figsize=(17, 15), gridspec_kw={'height_ratios': [0.5, 1.5, 1.5]})
    # fig.suptitle("Pressure-time conditions")

    # Plot 0: convergence rate
    ax[0].plot(cr["time"]/1e6, cr["conv_rate"], color="grey", linewidth=3, label = "Convergence rate (cm/yr)")
    ax[0].legend(loc='upper right', fontsize = labels_font)
    ax[0].tick_params(axis='both', which='major', labelsize=labels_font-2)
    # ax[0].set_ylabel("Convergence rate (cm/yr)", color="grey", fontweight="bold")
    ax[0].set_xlabel("Time (Ma)")
    ax[0].set_title("Convergence rate and exhumed particles")
    ax[0].set_ylim(0, cr["conv_rate"].max()+0.2)
    ax[0].set_xlim(0, 51)
    ax[0].label_outer()  # Only show outer labels and tick labels
    # vertical line at 35 Ma
    ax[0].axvline(x=35, color="grey", linestyle="--", label="35 Ma")

    # Plot 1: Pressure-time for exhumed particles
    legend_tmax = ax[1].legend(handles=handles_tmax, title="Peak Pressure", loc=(0.69, 0.02), frameon=True, fontsize = labels_font, title_fontsize=labels_font)
    ax[1].add_artist(legend_tmax)  # Manually add the tmax legend
    ax[1].legend(handles=handles_tfin, title="Time at Exhumability", loc=(0.814, 0.02), frameon=True, fontsize=labels_font, title_fontsize=labels_font)

    ax[1].add_patch(plt.Rectangle((0, 0), exhumed_list["tin"].max(), 2000, color='grey', alpha=0.5))
    ax[1].text(1, 0.5, " Particles \n Subduction", fontsize=18)
    sns.scatterplot(data=exhumed_list, x="tmax", y="maxP", hue="lithology", ax=ax[1], palette=colors_tmax, legend=False, s=80, linewidth=0.5, edgecolor="white")
    sns.scatterplot(data=exhumed_list, x="tfin", y="Pexh", hue="lithology", ax=ax[1], palette=colors_tfin, legend=False, s=80, linewidth=0.5, edgecolor="white")
    # ax[1].set_xlabel("Time (Ma)")
    ax[1].set_ylabel("Pressure (GPa)", fontsize=labels_font)
    #make ticklabels font bigger
    ax[1].tick_params(axis='both', which='major', labelsize=labels_font-2)
    ax[1].set_xlim(0, 51)
    ax[1].set_ylim(0, exhumed_list["maxP"].max()+0.1)
    ax[1].label_outer()  # Only show outer labels and tick labels
    ax[1].axvline(x=35, color="grey", linestyle="--", label="35 Ma")

    # Plot 3: Histogram of tmax and tfin
    ax[2].add_patch(plt.Rectangle((0, 0), exhumed_list["tin"].max(), 2000, color='grey', alpha=0.5))
    ax[2].text(1, 30, " Particles \n Subduction", fontsize=18)
    sns.histplot(
        data=exhumed_list, x="tmax", hue="lithology", ax=ax[2],
        palette=colors_tmax, legend=False, alpha=1, multiple="stack",
        edgecolor="black", linewidth=0.5
    )
    sns.histplot(
        data=exhumed_list, x="tfin", hue="lithology", ax=ax[2],
        palette=colors_tfin, legend=False, alpha=1, multiple="stack",
        edgecolor="black", linewidth=0.5
    )
    ax[2].set_yscale("log")
    ax[2].set_xlabel("Time (Ma)", fontsize=labels_font)
    ax[2].set_ylabel("Log(Particles count)", fontsize=labels_font)
    ax[2].tick_params(axis='both', which='major', labelsize=labels_font-2)
    ax[2].set_xlim(0, 51)
    ax[2].set_ylim(1, 2.e3)
    ax[2].axvline(x=35, color="grey", linestyle="--", label="35 Ma")



    # add no vertical space between the plots
    plt.subplots_adjust(hspace=0)
    

    plt.savefig(f"{plot_loc}/scatter_exhumVScr.png", dpi=500)
    plt.close()





    # # Scatter plots for exhumed particles
    # sns.scatterplot(data=exhumed_list, x="tin", y="Pin", hue="lithology", ax=ax[0], palette=colors_tin, legend=True)
    # sns.scatterplot(data=exhumed_list, x="tmax", y="maxP", hue="lithology", ax=ax[0], palette=colors_tmax, legend=False)
    # sns.scatterplot(data=exhumed_list, x="tfin", y="Pexh", hue="lithology", ax=ax[0], palette=colors_tfin, legend=False)
    # # sns.scatterplot(data=exhumed_list, x="tin", y="maxP", hue="lithology", ax=ax[0], palette=colors_tin, legend=True, size="vbur_category", sizes=sizes)
    # # sns.scatterplot(data=exhumed_list, x="tmax", y="maxP", hue="lithology", ax=ax[0], palette=colors_tmax, legend=False)
    # # sns.scatterplot(data=exhumed_list, x="tfin", y="maxP", hue="lithology", ax=ax[0], palette=colors_tfin, legend=False, style="vexh_category", markers=markers)
    # a1 = ax[0].twinx()
    # a1.plot(cr["time"]/1e6, cr["conv_rate"], color="grey", label="Convergence rate", linewidth=3)
    # a1.set_ylabel("Convergence rate (cm/yr)", color="grey", fontweight="bold")
    # ax[0].set_xlabel("Time (Ma)")
    # ax[0].set_ylabel("Pressure (GPa)")
    # ax[0].set_title("Exhumed particles")
    # ax[0].set_ylim(0, stagnant_list["maxP"].max()+0.2)

    # # Create a custom legend for marker sizes
    # handles, labels = ax[0].get_legend_handles_labels()
    # unique_labels = set(labels)
    # new_handles = []
    # new_labels = []

    # # Add lithology colors
    # # for label in unique_labels:
    # #     if label in colors_tin.keys():
    # #         new_handles.append(plt.Line2D([0], [0], marker='o', color='w', label=label,
    # #                                         markerfacecolor=colors_tin[label], markersize=10))
    # #         new_labels.append(label)

    # # # Add marker sizes for vbur categories with bin values
    # # for category, size in sizes.items():
    # #     low_value = e_bin_edges_bur[list(sizes.keys()).index(category)]
    # #     high_value = e_bin_edges_bur[list(sizes.keys()).index(category) + 1]
    # #     new_handles.append(plt.Line2D([0], [0], marker='o', color='w', label=f"vbur {category} ({low_value}, {high_value})",
    # #                                     markerfacecolor='black', markersize=size / 7))  # Scale down for consistency
    # #     new_labels.append(f"vbur {category} ({low_value}, {high_value})")

    # # # Add marker sizes for vexh categories with bin values
    # # for category, marker in markers.items():
    # #     low_value = e_bin_edges_exh[list(markers.keys()).index(category)]
    # #     high_value = e_bin_edges_exh[list(markers.keys()).index(category) + 1]
    # #     new_handles.append(plt.Line2D([0], [0], marker=marker, color='w', label=f"vexh {category} ({low_value}, {high_value})",
    # #                                     markerfacecolor='black', markersize=7))  # Scale down for consistency
    # #     new_labels.append(f"vexh {category} ({low_value}, {high_value})")

    # # ax[0].legend(handles=new_handles, labels=new_labels, title="Lithology and velocities", loc='upper right')



    # # Scatter plots for stagnant particles
    # # sns.scatterplot(data=stagnant_list, x="tin", y="maxP", hue="lithology", ax=ax[1], palette=colors_tin, legend=True, size="vbur_category", sizes=sizes)
    # # sns.scatterplot(data=stagnant_list, x="tmax", y="maxP", hue="lithology", ax=ax[1], palette=colors_tmax, legend=False)
    # # sns.scatterplot(data=stagnant_list, x="tfin", y="maxP", hue="lithology", ax=ax[1], palette=colors_tfin, legend=False, style="vstag_category", markers=markers)
    # sns.scatterplot(data=stagnant_list, x="tin", y="Pin", hue="lithology", ax=ax[1], palette=colors_tin, legend=True)
    # sns.scatterplot(data=stagnant_list, x="tmax", y="maxP", hue="lithology", ax=ax[1], palette=colors_tmax, legend=False)
    # sns.scatterplot(data=stagnant_list, x="tfin", y="Pstag", hue="lithology", ax=ax[1], palette=colors_tfin, legend=False)
    # a2 = ax[1].twinx()
    # a2.plot(cr["time"]/1e6, cr["conv_rate"], color="grey", label="Convergence rate", linewidth=3)
    # a2.set_ylabel("Convergence rate (cm/yr)", color="grey", fontweight="bold")
    # ax[1].set_xlabel("Time (Ma)")
    # ax[1].set_ylabel("Pressure (GPa)")
    # ax[1].set_title("Stagnant particles")
    # ax[1].set_ylim(0, stagnant_list["maxP"].max()+0.2)   
    # ax[1].legend(loc='upper center')

    # # Create a custom legend for marker sizes
    # handles, labels = ax[1].get_legend_handles_labels()
    # unique_labels = set(labels)
    # new_handles = []
    # new_labels = []

    # # # Add lithology colors
    # # for label in unique_labels:
    # #     if label in colors_tin.keys():
    # #         new_handles.append(plt.Line2D([0], [0], marker='o', color='w', label=label,
    # #                                         markerfacecolor=colors_tin[label], markersize=10))
    # #         new_labels.append(label)
    
    # # # Add marker sizes for vbur categories with bin values
    # # for category, size in sizes.items():
    # #     low_value = s_bin_edges_bur[list(sizes.keys()).index(category)]
    # #     high_value = s_bin_edges_bur[list(sizes.keys()).index(category) + 1]
    # #     new_handles.append(plt.Line2D([0], [0], marker='o', color='w', label=f"vbur {category} ({low_value}, {high_value})",
    # #                                     markerfacecolor='black', markersize=size / 7))  # Scale down for consistency
    # #     new_labels.append(f"vbur {category} ({low_value}, {high_value})")
    
    # # # Add marker sizes for vexh categories with bin values
    # # for category, marker in markers.items():
    # #     low_value = s_bin_edges_stag[list(markers.keys()).index(category)]
    # #     high_value = s_bin_edges_stag[list(markers.keys()).index(category) + 1]
    # #     new_handles.append(plt.Line2D([0], [0], marker=marker, color='w', label=f"vstag {category} ({low_value}, {high_value})",
    # #                                     markerfacecolor='black', markersize=7))  # Scale down for consistency
    # #     new_labels.append(f"vstag {category} ({low_value}, {high_value})")
    
    # # ax[1].legend(handles=new_handles, labels=new_labels, title="Lithology and velocities", loc=(0.5, 0.37))

    # format = ['png', 'eps']
    # for f in format:
    #     plt.savefig(f"{plot_loc}/pressure_time.{f}", dpi=500)
    # # plt.savefig(f"{plot_loc}/pressure_time.{format}", dpi=500)
    # plt.close()

if __name__ == "__main__":
    main()

