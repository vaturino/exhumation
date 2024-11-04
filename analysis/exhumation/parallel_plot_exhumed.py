#! /usr/bin/python3

### Script that loops through all of the particles and identifies and plots the exhumed particles ###

import pandas as pd
import numpy as np
import json
from tqdm import tqdm
import matplotlib.pyplot as plt
import argparse
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import seaborn as sns
from matplotlib.cm import get_cmap

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
pd.options.mode.chained_assignment = None  # default='warn'

path = str(Path(Path(__file__).parent.absolute()).parent.absolute())
sys.path.insert(0, path)
from libraries.functions import *
from libraries.particles import *
from libraries.exhumation import *  


import numpy as np

def process_particle(p, pt_files, line_colors, compositions, composition_mapping, ymax = 900.):
    pt_single = pd.read_csv(f"{pt_files}/pt_part_{p}.txt", sep="\s+")

    # pt_single["Plith"] = (ymax - pt_single["depth"]) * 1e3 * 9.81 * 3100 / 1e9
    pt_single["T"] = (pt_single["T"]) # + 0.6*(ymax - pt_single["depth"])) - 273.
    pt_single["terrain"] = sum((ind_c + 1) * pt_single[c].round() for ind_c, c in enumerate(compositions))
    pt_single["terrain"] = pt_single["terrain"].round()
    pt_single["lithology"] = pt_single["terrain"].map(composition_mapping)
    pt_single["lithology"] = pt_single["lithology"].iloc[-1]

    max_depth = ymax - pt_single["depth"].min()
   
    stagnant = ((pt_single.depth.iat[-1] - pt_single.depth.min()) > 0.01) and ((pt_single.depth.iat[-1] - pt_single.depth.min()) < 0.25 * max_depth) and (max_depth > 10.)
    exhumed = ((pt_single.depth.iat[-1] - pt_single.depth.min()) >= 0.25 * max_depth) and (max_depth > 10.)

    result = {
        "p": p,
        "subducted": pt_single["P"].max() > 3.,
        "exhumed": exhumed, 
        "stagnant": stagnant,
        "pt_single": pt_single,
        "max_depth": max_depth,
        "line_color": line_colors[p]
    }

    if not result["subducted"] and not result["exhumed"] and not result["stagnant"]:
        result["subducted"] = True

    return result



def main(): 
    parser = argparse.ArgumentParser(description='Script that gets some models and gives the kinematic indicators')
    parser.add_argument('json_file', help='json file with model name, time at the end of computation, folder where to save the plots, legend')
    args = parser.parse_args()

    models_loc = '/home/vturino/Vale_nas/exhumation/raw_outputs/'
    csvs_loc = '/home/vturino/Vale_nas/exhumation/gz_outputs/'
    json_loc = '/home/vturino/PhD/projects/exhumation/pyInput/'

    with open(f"{json_loc}{args.json_file}") as json_file:
        configs = json.load(json_file)
    compositions = configs['compositions']

    plot_loc = f"/home/vturino/PhD/projects/exhumation/plots/single_models/{configs['models'][0]}"
    txt_loc = f'{plot_loc}/txt_files'
    os.makedirs(txt_loc, exist_ok=True)

    eloc = f"{plot_loc}/exhumed"
    os.makedirs(eloc, exist_ok=True)

    exhloc = f"{txt_loc}/exhumed"
    os.makedirs(exhloc, exist_ok=True)

    pt_files = f'{txt_loc}/PT'
    npa = len(os.listdir(pt_files))
    print("Total number of particles = ", npa)

    # trench_pos = pd.read_csv(f"{txt_loc}/trench_pos.txt", sep="\s+")

    init = pd.read_csv(f"{txt_loc}/particles_indexes.txt", sep="\s+")
    pal1 = plt.get_cmap('viridis')
    norm = plt.Normalize(init["init_x"].min() / 1e3, init["init_x"].max() / 1e3)
    line_colors = pal1(norm(init["init_x"] / 1e3))

    exh = pd.DataFrame(columns=["id", "maxPP", "maxPT", "maxTP", "maxTT", "lithology", "tmax"], index=range(npa))
    stag = pd.DataFrame(columns=["id", "maxPP", "maxPT", "maxTP", "maxTT", "lithology", "tmax"], index=range(npa))
    subd = pd.DataFrame(columns=["id"], index=range(npa))
    s_m = plt.cm.ScalarMappable(cmap=pal1, norm=norm)
    s_m.set_array([])

    subducted = 0
    exhumed = 0
    stagnant = 0
    

    composition_mapping = {ind_c + 1: c for ind_c, c in enumerate(compositions)}

    f1, a1 = plt.subplots(1, 2, figsize=(15, 5))
    f2, a2 = plt.subplots(1, 1)
    f3, a3 = plt.subplots(1, 2, figsize=(15, 5))

    with ProcessPoolExecutor() as executor:
        futures = [executor.submit(process_particle, p, pt_files, line_colors, compositions, composition_mapping) for p in range(npa)]
        for future in tqdm(as_completed(futures), total=npa):
            result = future.result()
            p = result["p"]
            pt_single = result["pt_single"]
            line_color = result["line_color"]

            ymax = 900.

            if result["subducted"]:
                subducted += 1
                a2.plot(pt_single["T"] , pt_single["Plith"], color=line_color)
                subd.iloc[p] = [p]
            elif result["exhumed"]:
                exhumed += 1
                exh.loc[p] = [p, pt_single["Plith"].max(), pt_single["T"].iloc[pt_single["Plith"].idxmax()]  , pt_single["Plith"].iloc[pt_single["T"].idxmax()], pt_single["T"].iloc[pt_single["T"].idxmax()] , pt_single["lithology"].iloc[-1], pt_single["time"].iloc[pt_single["Plith"].idxmax()] / 2]
                a1[0].plot(pt_single["T"] , pt_single["Plith"], color=line_color)
                a1[1].plot((pt_single["x"] )/ 1.e3, ymax - pt_single["depth"], color=line_color)
                a1[1].invert_yaxis()
            elif result["stagnant"]:
                stagnant += 1
                stag.loc[p] = [p, pt_single["Plith"].max(), pt_single["T"].iloc[pt_single["Plith"].idxmax()] , pt_single["Plith"].iloc[pt_single["T"].idxmax()], pt_single["T"].iloc[pt_single["T"].idxmax()] , pt_single["lithology"].iloc[-1], pt_single["time"].iloc[pt_single["Plith"].idxmax()] / 2]
                a3[0].plot(pt_single["T"] , pt_single["Plith"], color=line_color)
                a3[1].plot(pt_single["x"] / 1.e3, ymax - pt_single["depth"], color=line_color)
                a3[1].set_ylim(72,-2)
                # a3[1].invert_yaxis()

    for ax, title in zip([a1[0], a1[1], a2, a3[0], a3[1]], ["Exhumed particles", "Exhumed particles trajectory", "Subducted particles", "Stagnant particles", "Stagnant particles trajectory"]):
        ax.set_title(title)
    a1[0].set_xlabel("Temperature (C)")
    a1[0].set_ylabel("Pressure (GPa)")
    a1[1].set_xlabel("Distance (km)")
    a1[1].set_ylabel("Depth (km)")
    a2.set_xlabel("Temperature (C)")
    a2.set_ylabel("Pressure (GPa)")
    a2.set_ylim(0, 4.5)
    a2.set_xlim(0, 1100)
    a3[0].set_xlabel("Temperature (C)")
    a3[0].set_ylabel("Pressure (GPa)")
    a3[1].set_xlabel("Distance (km)")
    a3[1].set_ylabel("Depth (km)")

    for fig, ax in zip([f1, f2, f3], [a1[1], a2, a3[1]]):
        plt.colorbar(s_m, ax=ax).set_label('Initial X Position (km)')
        fig.tight_layout()

    f1.savefig(f"{plot_loc}/possibly_exhumed.png")
    f3.savefig(f"{plot_loc}/stagnant.png")
    f2.savefig(f"{plot_loc}/filtered_out.png", dpi = 500)
    plt.close()

    exh.dropna(inplace=True)
    stag.dropna(inplace=True)
    print("Subducted particles = ", subducted)
    print("Exhumed particles = ", len(exh))
    print("Stagnant particles = ", len(stag))

     #Titles 
    exh_title_str = ""
    for j in range(len(exh.lithology.value_counts())):
        lithology = exh.lithology.value_counts().index[j]
        count = exh.lithology.value_counts().values[j]
        percentage = (count / npa) * 100
        exh_title_str += f"{lithology}: count = {count}, percentage = {percentage:.5f}%\n"

    stag_title_str = ""
    for j in range(len(stag.lithology.value_counts())):
        lithology = stag.lithology.value_counts().index[j]
        count = stag.lithology.value_counts().values[j]
        percentage = (count / npa) * 100
        stag_title_str += f"{lithology}: count = {count}, percentage = {percentage:.5f}%\n"


    def plot_max_conditions(df, title_str, plot_file, lith_count):
        f, a = plt.subplots(2, 2, figsize=(10, 10))

        # Get unique lithologies and assign a color palette
        unique_lithologies = df.lithology.unique()
        colors = sns.color_palette("colorblind", len(unique_lithologies))

        # Create a mapping of lithologies to colors
        color_mapping = {lithology: color for lithology, color in zip(unique_lithologies, colors)}

        # Scatterplot with consistent colors
        sns.scatterplot(ax=a[0, 0], data=df, x="maxPT", y="maxPP", hue="lithology", palette=color_mapping, linewidth=0.2, size="tmax")
        a[0, 0].set_xlabel("T ($^\circ$C)")
        a[0, 0].set_ylabel("P (GPa)")
        a[0, 0].set_title("Max P")

        sns.scatterplot(ax=a[0, 1], data=df, x="maxTT", y="maxTP", hue="lithology", palette=color_mapping, linewidth=0.2, size="tmax")
        a[0, 1].set_xlabel("T ($^\circ$C)")
        a[0, 1].set_ylabel("P (GPa)")
        a[0, 1].set_title("Max T")

        # Pie chart with consistent colors
        lithology_counts = df.lithology.value_counts()
        a[1, 0].pie(lithology_counts, labels=lithology_counts.index, autopct='%1.1f%%', 
                    colors=[color_mapping[lithology] for lithology in lithology_counts.index])
        a[1, 0].set_title("Lithology")

        # Histogram with hue based on lithology
        sns.histplot(ax=a[1, 1], data=df, x="tmax", bins=20, hue="lithology", element="step", palette=color_mapping)
        a[1, 1].set_xlabel("Peak pressure Time (Ma)")

        f.suptitle(title_str)
        f.text(0.2, 0.05, lith_count, ha='center', va='center')
        f.tight_layout()
        plt.savefig(plot_file, dpi=1000)
        plt.close()


    plot_max_conditions(exh, f"Total number of exhumed particles = {len(exh)} - {(len(exh) / npa) * 100:.1f}%", f"{plot_loc}/max_PT_conditions.png", exh_title_str)
    plot_max_conditions(stag, f"Total number of stagnant particles = {len(stag)} - {(len(stag) / npa) * 100:.1f}%", f"{plot_loc}/max_PT_conditions_stagnant.png", stag_title_str)

    exh.to_csv(f"{txt_loc}/exhumed_particles.txt", sep="\t", index=False)
    stag.to_csv(f"{txt_loc}/stagnant_particles.txt", sep="\t", index=False)
    subd.to_csv(f"{txt_loc}/subducted_particles.txt", sep="\t", index=False)

if __name__ == "__main__":
    main()
