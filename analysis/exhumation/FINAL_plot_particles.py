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


# load the pressure-temperature data from the txt files extracted through the "parallel_get_PT.py" script
def load_data(id, txt_loc):
    data = pd.read_csv(f"{txt_loc}/PT/pt_part_{id}.txt", sep="\s+")
    data["time"] = data["time"] / 2
    data["ts"] = data["time"] * 2
    return data


# Compute the time intervals for the particles
def compute_time_intervals(data, stagnation_min, time_thresh):
    start_time = None
    data["time_interval"] = np.nan 
    for i, row in data.iterrows():
        if row['duration'] == time_thresh:
            if start_time is None:
                start_time = row['time'] 
            end_time = row['time']

            if i == len(data) - 1 or data.at[i+1, 'duration'] != time_thresh:
                for j in range(i, -1, -1):
                    if data.at[j, 'duration'] == time_thresh:
                        if start_time != end_time:  
                            data.at[j, 'time_bin'] = f"[{start_time}, {end_time})"
                            data.at[j, 'ti'] = start_time
                            data.at[j, 'tf'] = end_time
                            data.at[j, 'time_interval'] = end_time - start_time
                    else:
                        break
                start_time = None  
            else:
                continue
        else:
            start_time = None

    data.loc[data["time_interval"] < stagnation_min, ["ti", "tf", "time_bin", "time_interval"]] = np.nan

    return data


def calculate_middle_values(data):
    avg_values = data.groupby('time_bin').agg({'Plith': 'mean', 'time': 'mean', 'T': 'mean'}).reset_index()
    data["Pm"] = np.nan
    data["tm"] = np.nan
    data["Tm"] = np.nan

    for tint in avg_values["time_bin"]:
        data["Pm"] = np.where(data["time_bin"] == tint, avg_values[avg_values["time_bin"] == tint]["Plith"].values[0], data["Pm"])
        data["tm"] = np.where(data["time_bin"] == tint, avg_values[avg_values["time_bin"] == tint]["time"].values[0], data["tm"])
        data["Tm"] = np.where(data["time_bin"] == tint, avg_values[avg_values["time_bin"] == tint]["T"].values[0], data["Tm"])

    return data
    
def process_times_for_particle(data, stagnation_min, time_thresh, grad_thresh):
    data['gradient'] = np.gradient(data["Plith"], data["time"])
    lowgrad = data[abs(data["gradient"]) < grad_thresh].reset_index(drop=True)
    lowgrad["duration"] = lowgrad["time"].diff()
    lowgrad['time_bin'] = None


    lowgrad = compute_time_intervals(lowgrad, stagnation_min, time_thresh)
    lowgrad = calculate_middle_values(lowgrad)
    lowgrad = lowgrad[lowgrad["Plith"] < 6.0]
    return lowgrad


def assign_particle_values(particles, lowgrad_dyn, lowgrad_kin, lowgrad_trans, i, c):
            for col in c:
                if not lowgrad_dyn[col].empty:
                    particles[f"{col}_dyn"].iloc[i] = lowgrad_dyn[col].values[0]
                else:
                    particles[f"{col}_dyn"].iloc[i] = np.nan

                if not lowgrad_kin[col].empty:
                    particles[f"{col}_kin"].iloc[i] = lowgrad_kin[col].values[0]
                else:
                    particles[f"{col}_kin"].iloc[i] = np.nan

                if not lowgrad_trans[col].empty:
                    particles[f"{col}_trans"].iloc[i] = lowgrad_trans[col].values[0]
                else:
                    particles[f"{col}_trans"].iloc[i] = np.nan

            return particles



def process_particle(p, txt_loc, line_colors, compositions, composition_mapping, thresh, time_thresh, stagnation_min, c, grad_thresh, ymax=901.):
    # Load data
    pt_single = load_data(p, txt_loc)

    # Process terrain and lithology
    pt_single["T"] = pt_single["T"]
    pt_single["terrain"] = sum((ind_c + 1) * pt_single[c].round() for ind_c, c in enumerate(compositions))
    pt_single["terrain"] = pt_single["terrain"].round()
    pt_single["lithology"] = pt_single["terrain"].map(composition_mapping)
    # pt_single["lithology"] = pt_single["lithology"].iloc[-1]
    

    max_depth = ymax - pt_single["depth"].min()
    exhumed = ((pt_single.depth.iat[-1] - pt_single.depth.min()) >= thresh * max_depth) and (max_depth > 10.)

    stagnant = not exhumed # Default as stagnant
    subducted = not exhumed and not stagnant  # Default as not subducted

    particle_data = None

    if stagnant:
        pt_single["Plith"] = pt_single["Plith"].rolling(window=10, min_periods=1).mean()
        pt_single["time"] = pt_single["time"].rolling(window=10, min_periods=1).mean()
        pt_single["gradient"] = np.gradient(pt_single["Plith"], pt_single["time"])
        lowgrad = process_times_for_particle(pt_single, stagnation_min, time_thresh, grad_thresh)

        bin_number = lowgrad["time_bin"].nunique()
        if bin_number == 0:
            subducted = True  # Mark particle as subducted if no bins are found
        else:
            # Process bins only if they exist
            lowgrad_dyn = lowgrad[lowgrad["tm"] < 33.]
            lowgrad_kin = lowgrad[lowgrad["tm"] > 37.]
            lowgrad_trans = lowgrad[(lowgrad["tm"] >= 33.) & (lowgrad["tm"] <= 37.)]

            particle_data = {
                "id": p,
                "lithology": pt_single["lithology"]
            }

            for col in c:
                particle_data[f"{col}_dyn"] = lowgrad_dyn[col].values[0] if not lowgrad_dyn[col].empty else np.nan
                particle_data[f"{col}_kin"] = lowgrad_kin[col].values[0] if not lowgrad_kin[col].empty else np.nan
                particle_data[f"{col}_trans"] = lowgrad_trans[col].values[0] if not lowgrad_trans[col].empty else np.nan
    

    # Determine the color value safely
    color = line_colors[p] if isinstance(line_colors, dict) else line_colors[p] if len(line_colors) > p else None

    result = {
        "p": p,
        "subducted": subducted,
        "exhumed": exhumed,
        "stagnant": stagnant,
        "pt_single": pt_single,
        "max_depth": max_depth,
        "line_color": color,
        "particle_data": particle_data,
        "data": particle_data if stagnant else None
    }

    return result




# Plot the maximum conditions for each particle
def plot_max_conditions(df, title_str, plot_file, lith_count, P, T, t):
        print(f"Plotting exhumed")
        f, a = plt.subplots(1, 3, figsize=(15, 5))

        
        for p in range(len(P)):
            # Scatterplot with consistent colors
            # Get unique lithologies and assign a color palette
            # print(df[lith[p]].unique())
            unique_lithologies = df.lithology.unique()
            colors = sns.color_palette("colorblind", len(unique_lithologies))
    
            # Create a mapping of lithologies to colors
            color_mapping = {lithology: color for lithology, color in zip(unique_lithologies, colors)}

            sns.scatterplot(ax=a[0], data=df, x=T[p], y=P[p], hue="lithology", palette=color_mapping, linewidth=0.2)
            a[0].set_xlabel("T ($^\circ$C)")
            a[0].set_ylabel("P (GPa)")
            a[0].set_title("Max P")


            # Pie chart with consistent colors
            lithology_counts = df.lithology.value_counts()
            a[1].pie(lithology_counts, labels=lithology_counts.index, autopct='%1.1f%%', 
                        colors=[color_mapping[lithology] for lithology in lithology_counts.index])
            a[1].set_title("Lithology")

            # Histogram with hue based on lithology
            sns.histplot(ax=a[2], data=df, x=t[p], bins=20, hue="lithology", element="step", palette=color_mapping)
            a[2].set_xlabel("Peak pressure Time (Ma)")

            f.suptitle(title_str)
            f.text(0.5, 0.1, lith_count, ha='center', va='center')
            f.tight_layout()
        plt.savefig(plot_file, dpi=1000)
        plt.close()





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

    stagnation_min = 10.
    grad_thresh = 0.01
    time_thresh = 0.5
    exhumed_thresh = 0.25
    c = ["Pm", "tm", "Tm", "time_interval", "ti", "tf", "lithology"]


    init = pd.read_csv(f"{txt_loc}/particles_indexes.csv")
    pal1 = plt.get_cmap('viridis')
    norm = plt.Normalize(init["init_x"].min() / 1e3, init["init_x"].max() / 1e3)
    line_colors = pal1(norm(init["init_x"] / 1e3))

    exh = pd.DataFrame(columns=["id", "maxPP", "maxPT", "maxTP", "maxTT", "lithology", "tmax"], index=range(npa))
    subd = pd.DataFrame(columns=["id"], index=range(npa))

    fixed_columns = ["id", "lithology"]
    timing = ["dyn", "trans", "kin"]
    columns = fixed_columns + [f"{c}_{t}" for c in ["Pm", "tm", "Tm", "time_interval", "ti","tf", "lithology"] for t in timing]
    stag = pd.DataFrame(columns=columns, index=range(npa))


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
        futures = [
            executor.submit(
                process_particle, 
                p, txt_loc, line_colors, compositions, composition_mapping, exhumed_thresh, time_thresh, stagnation_min, c, grad_thresh
            ) 
            for p in range(npa)
        ]
        for future in tqdm(as_completed(futures), total=npa):
            try:
                result = future.result()
                if result is None:
                    continue  # Skip to the next future if the result is None
                p = result["p"]
                pt_single = result["pt_single"]
                line_color = result["line_color"]

                if result["subducted"]:
                    subducted += 1
                    a2.plot(pt_single["T"], pt_single["Plith"], color=line_color)
                    subd.iloc[p] = [p]
                elif result["exhumed"]:
                    exhumed += 1
                    exh.loc[p] = [
                        p, 
                        pt_single["Plith"].max(), 
                        pt_single["T"].iloc[pt_single["Plith"].idxmax()],  
                        pt_single["Plith"].iloc[pt_single["T"].idxmax()],
                        pt_single["T"].iloc[pt_single["T"].idxmax()],
                        pt_single["lithology"].iloc[-1],
                        pt_single["time"].iloc[pt_single["Plith"].idxmax()]
                    ]
                    a1[0].plot(pt_single["T"], pt_single["Plith"], color=line_color)
                    a1[1].plot((pt_single["time"]), pt_single["Plith"], color=line_color)
                    a1[1].set_xlim(0, 50)
                elif result["stagnant"]:
                    stagnant += 1
                    if result["data"] is not None:
                        result_data_cleaned = result["data"].copy()
                        result_data_cleaned["lithology"] = result_data_cleaned["lithology"].iloc[-1] 
                        particle_data_series = pd.Series(result_data_cleaned, index=stag.columns)
                        stag.loc[p] = particle_data_series          
                    a3[0].plot(pt_single["T"], pt_single["Plith"], color=line_color)
                    a3[1].plot((pt_single["time"]), pt_single["Plith"], color=line_color)
                    a3[1].set_xlim(0, 50)

            except Exception as e:
                print(f"Error processing particle {p}: {e}")




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
    columns_to_check = [f"{col}_dyn" for col in c] + [f"{col}_kin" for col in c] + [f"{col}_trans" for col in c]
    stag.dropna(subset=columns_to_check, how='all', inplace=True)
    subd.dropna(inplace=True)

    print("Subducted particles = ", len(subd))
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


    # Calculate the highest numerical value between Pm_kin, Pm_trans, Pm_dyn
    stag["Pm"] = stag[["Pm_kin", "Pm_trans", "Pm_dyn"]].max(axis=1)
    # max_Pm_column = stag[["Pm_kin", "Pm_trans", "Pm_dyn"]].idxmax(axis=1)
    # # Map each column to its corresponding tm_i column
    # tm_mapping = {
    #     "Pm_kin": "tm_kin",
    #     "Pm_trans": "tm_trans",
    #     "Pm_dyn": "tm_dyn"
    # }


    

    plot_max_conditions(exh, f"Total number of exhumed particles = {len(exh)} - {(len(exh) / npa) * 100:.1f}%", f"{plot_loc}/max_PT_conditions.png", exh_title_str, ["maxPP"], ["maxPT"], ["tmax"])
    plot_max_conditions(stag, f"Total number of stagnant particles = {len(stag)} - {(len(stag) / npa) * 100:.1f}%", f"{plot_loc}/max_PT_conditions_stagnant.png", stag_title_str, ["Pm_kin", "Pm_dyn", "Pm_trans"], ["Tm_kin", "Tm_dyn", "Tm_trans"], ["tm_kin", "tm_dyn", "tm_trans"])   
    print("particles plotted")

    exh.to_csv(f"{txt_loc}/exhumed_particles.txt", sep="\t", index=False)
    stag = stag.astype({col: 'float' for col in columns if col not in fixed_columns+["lithology_kin", "lithology_dyn", "lithology_trans"]})
    stag.to_csv(f"{txt_loc}/stagnant_particles.txt", sep=" ", index=False, float_format='%.2f', na_rep="NaN")
    subd.to_csv(f"{txt_loc}/subducted_particles.txt", sep="\t", index=False)

    print("All particles processed and saved")

if __name__ == "__main__":
    main()
