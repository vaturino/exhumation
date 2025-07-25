#!/usr/bin/env python3

import pandas as pd
import numpy as np
import json
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import os

import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42  # TrueType fonts, keeps text as text
mpl.rcParams['ps.fonttype'] = 42   # In case you ever save EPS
mpl.rcParams['text.usetex'] = False  # Avoid text being converted to paths



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

    # Colors
    color = ["mediumslateblue", "firebrick", "darkslateblue", "forestgreen", "darkgreen"]
    colors_tfin = {
        "sed": "cornflowerblue",
        "oc": "darkorange",
        "ecl": "#61aa2f",
        "serp": "lightsalmon"
    }

    fig, ax = plt.subplots(len(configs['models']) + 1, 1, figsize=(8, 12), sharex=True, height_ratios=[1] * len(configs['models']) + [0.35])

    for m, mod in enumerate(configs['models']):
        mod_loc = f"/home/vturino/PhD/projects/exhumation/plots/single_models/{mod}"
        exhumed_list = pd.read_csv(f"{mod_loc}/txt_files/exhumed_times.txt", sep="\s+")

        bin_edges = np.arange(0, 51, 1)  # Define bin edges

        # Sort the exhumed_list DataFrame to ensure "oc" lithology is plotted in front
        exhumed_list['lithology'] = pd.Categorical(
            exhumed_list['lithology'], 
            categories=["oc", "sed", "ecl", "serp"], 
            ordered=True
        )
        exhumed_list = exhumed_list.sort_values(by="lithology")

        # Plot exhumed particles with hue as lithology and color by tfin
        sns.histplot(
            exhumed_list,
            x="texh",
            hue="lithology",
            ax=ax[m],
            bins=bin_edges,
            palette=colors_tfin,
            alpha=0.5,
            linewidth=1.3,
            element="step",
            legend=False,
        )
        ax[m].axvline(x=35, color="grey", linestyle="--", linewidth=1)
        ax[m].set_yscale("log")
        ax[m].set_xlim(0, 50)
        ax[m].set_ylim(0.7, 5.e3)
        ax[m].set_ylabel("Count (log)")
        ax[m].text(0.05, 0.8, configs["names"][m], horizontalalignment='left', verticalalignment='center', transform=ax[m].transAxes)

    # Add the convergence rate vs time plot for the reference model
    ref_model = configs['models'][1]
    ref_loc = f"/home/vturino/PhD/projects/exhumation/plots/single_models/{ref_model}"
    conv_rate_data = pd.read_csv(f"{ref_loc}/txt_files/2D_v.txt", sep="\s+")
    conv_rate_data.iloc[0] = conv_rate_data.iloc[1]  # Replace the first value with the second to avoid NaN

    ax[-1].plot(conv_rate_data['time']/1.e6, conv_rate_data['conv_rate'], color="k")
    ax[-1].axvline(x=35, color="grey", linestyle="--", linewidth=1)
    ax[-1].set_xlim(0, 50)
    ax[-1].set_ylim(0, 8)
    ax[-1].set_ylabel("$v_{c} (cm/yr)$")
    ax[-1].set_xlabel("Time (Ma)")
    # ax[-1].legend(loc="upper right")

    # ax[0].legend(loc='upper right', title="Lithology")

    plt.suptitle("Time at the exhumation threshold", fontsize=16)

    plt.tight_layout()
    plt.savefig(f"{plotloc}/rheo_exhumed_particles_with_conv_rate.pdf")
    plt.close()


if __name__ == "__main__":
    main()
