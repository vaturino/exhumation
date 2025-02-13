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
    exhumed_list = pd.read_csv(f"{txt_loc}/exhumed_particles.txt", sep="\s+")
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
    fig,ax = plt.subplots(1, 2, figsize=(10, 5))



    # # Histograms for tmax and tfin (Manual layering with plt.bar)
    bin_edges = np.arange(0, 50,1)  # Define bin edges


    # Sort exhumed_list by lithology to ensure 'oc' is plotted over 'sed'
    exhumed_list_sorted = exhumed_list.sort_values(by='lithology', ascending=True)
    sns.histplot(exhumed_list_sorted, x="tmax", hue="lithology", bins=bin_edges, ax=ax[0], palette=colors_tfin, alpha=alpha, linewidth=1, element="step", legend=False)
    

    # Set scale, labels, and limits
    ax[0].set_yscale("log")
    ax[0].set_xlabel("Time (Myr)")
    ax[0].set_ylabel("Particles count (log)")
    ax[0].set_xlim(0, 50)
    ax[0].set_ylim(0.7, 5.e3)
    ax[0].set_title("Exhumed particles: time at peak P")


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


    sns.histplot(stagnant_list_expanded, x="ti", ax=ax[1], hue = "lith_time", palette= colors_tfin, alpha=alpha, linewidth=1, element="step", bins= bin_edges, legend=False)
    
    
    ax[1].text(17, 2400, "Beginning of stagnation", fontsize=12)

    # Add a title, labels, and legend
    ax[1].set_xlabel("Time (Myr)")
    ax[1].set_ylabel("")
    ax[1].set_xlim(0, 50)
    ax[1].set_ylim(0.7, 5.e3)
    ax[1].set_yscale("log")
    ax[1].set_title("Stagnant particles: time at beginning of stagnation")






   
    plt.tight_layout()





    

    plt.savefig(f"{plot_loc}/peak_vs_stag.pdf", dpi=500)
    plt.close()





if __name__ == "__main__":
    main()
