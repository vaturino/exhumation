#! /usr/bin/python3

from matplotlib.cm import get_cmap
import pandas as pd
import numpy as np
import json
from tqdm import tqdm
import matplotlib.cm as cm
import matplotlib as mpl
import sys, os, subprocess
from matplotlib.pyplot import copper, show, xlabel, ylabel
import matplotlib.pyplot as plt
import argparse
from matplotlib.gridspec import GridSpec
from matplotlib import colors
import math as math
from scipy.signal import savgol_filter 
import seaborn as sns
from matplotlib.axes._axes import _log as matplotlib_axes_logger
matplotlib_axes_logger.setLevel('ERROR')
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
pd.options.mode.chained_assignment = None  # default='warn'
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

path = str(Path(Path(__file__).parent.absolute()).parent.absolute())
sys.path.insert(0, path)
from libraries.functions import *
from libraries.particles import *
from libraries.exhumation import *  
import plotly.express as px


############### FUNCTIONS ####################

def load_data(id, txt_loc):
    data = pd.read_csv(f"{txt_loc}/PT/pt_part_{id}.txt", sep="\s+")
    data["time"] = data["time"] / 2
    data["ts"] = data["time"] * 2
    return data


def compute_time_intervals(data, stagnation_min, time_thresh):
    # Determine the maximum pressure and its corresponding time
    Pmax = data["Plith"].max()
    tmax = data[data["Plith"] == Pmax]["time"].values[0]
    Pthresh = round(0.75 * Pmax, 3)
    
    # Find the exhumation threshold time (texh)
    texh = data[(data["Plith"] <= Pthresh) & (data["time"] > tmax)]["time"].min()
    if pd.isna(texh):  # Handle cases where texh cannot be determined
        texh = tmax - 0.5

    start_time = None  # Initialize start time for intervals

    # Iterate over the data rows
    for i, row in data.iterrows():
        if row['duration'] == time_thresh:  # Check if the row duration matches the threshold
            if start_time is None:
                start_time = row['time']
            end_time = row['time']

            # If the interval ends here, process it
            if i == len(data) - 1 or data.at[i + 1, 'duration'] != time_thresh and end_time ==50.:
                # Iterate backwards to assign interval values to relevant rows
                for j in range(i, -1, -1):
                    if data.at[j, 'duration'] == time_thresh:
                        # Adjust interval to start at texh if it spans across it
                        adjusted_start_time = max(start_time, texh)
                        
                        # Assign interval only if it begins at or after texh
                        if adjusted_start_time < end_time:
                            data.at[j, 'time_bin'] = f"[{adjusted_start_time}, {end_time})"
                            data.at[j, 'ti'] = adjusted_start_time
                            data.at[j, 'tf'] = end_time
                            data.at[j, 'time_interval'] = end_time - adjusted_start_time
                    else:
                        break  # Stop if the current row does not match the duration condition
                start_time = None  # Reset start time for the next interval
        else:
            start_time = None  # Reset start time if duration condition is not met

    return data



# def compute_time_intervals(data, stagnation_min, time_thresh):
#     # Determine the maximum pressure and its corresponding time
#     Pmax = data["Plith"].max()
#     tmax = data[data["Plith"] == Pmax]["time"].values[0]
#     Pthresh = round(0.75 * Pmax, 3)
    
#     # Find the exhumation threshold time (texh)
#     texh = data[(data["Plith"] <= Pthresh) & (data["time"] > tmax)]["time"].min()
#     if pd.isna(texh):  # Handle cases where texh cannot be determined
#         texh = tmax - 0.5

#     start_time = None  # Initialize start time for intervals

#     # Iterate over the data rows
#     for i, row in data.iterrows():
#         if row['duration'] == time_thresh:  # Check if the row duration matches the threshold
#             if start_time is None:
#                 start_time = row['time']
#             end_time = row['time']

#             # If the interval ends here, process it
#             if i == len(data) - 1 or data.at[i + 1, 'duration'] != time_thresh and end_time ==50.:
#                 # Iterate backwards to assign interval values to relevant rows
#                 for j in range(i, -1, -1):
#                     if data.at[j, 'duration'] == time_thresh:
#                         # Adjust interval to start at texh if it spans across it
#                         adjusted_start_time = max(start_time, texh)
                        
#                         # Assign interval only if it begins at or after texh
#                         if adjusted_start_time < end_time:
#                             data.at[j, 'time_bin'] = f"[{adjusted_start_time}, {end_time})"
#                             data.at[j, 'ti'] = adjusted_start_time
#                             data.at[j, 'tf'] = end_time
#                             data.at[j, 'time_interval'] = end_time - adjusted_start_time
#                     else:
#                         break  # Stop if the current row does not match the duration condition
#                 start_time = None  # Reset start time for the next interval
#         else:
#             start_time = None  # Reset start time if duration condition is not met


#     return data




def calculate_middle_values(data):
    # Calculate the average values for each time_bin
    avg_values = data.groupby('time_bin').agg({'Plith': 'mean', 'time': 'mean', 'T': 'mean'}).reset_index()

    # Initialize the columns
    data["Pm"] = np.nan
    data["tm"] = np.nan
    data["Tm"] = np.nan

    # Loop through each time_bin and assign the average values
    for tint in avg_values["time_bin"]:
        avg_Pm = avg_values[avg_values["time_bin"] == tint]["Plith"].values[0]
        avg_tm = avg_values[avg_values["time_bin"] == tint]["time"].values[0]
        avg_Tm = avg_values[avg_values["time_bin"] == tint]["T"].values[0]

        # Ensure these values are assigned to the truncated time intervals
        data.loc[data["time_bin"] == tint, "Pm"] = avg_Pm
        data.loc[data["time_bin"] == tint, "tm"] = avg_tm
        data.loc[data["time_bin"] == tint, "Tm"] = avg_Tm

    return data

    
def process_times_for_particle(data, stagnation_min, time_thresh, grad_thresh):
    data['gradient'] = np.gradient(data["Plith"], data["time"])
    lowgrad = data[abs(data["gradient"]) < grad_thresh].reset_index(drop=True)
    lowgrad["duration"] = lowgrad["time"].diff()
    lowgrad['time_bin'] = None

    lowgrad = compute_time_intervals(lowgrad, stagnation_min, time_thresh)
    lowgrad = calculate_middle_values(lowgrad)
    return lowgrad


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


def process_particle(id, part, txt_loc, stagnation_min, time_thresh, sloc, sfiles, grad_thresh):
    data = load_data(id, txt_loc)
    data['gradient'] = np.gradient(data["Plith"], data["time"])

    lowgrad = process_times_for_particle(data, stagnation_min, time_thresh, grad_thresh)
    lowgrad["lithology"] = part["lithology"]

    bin_number = lowgrad["time_bin"].nunique()
    if bin_number == 0:
        return None

    particle_data = {
        "id": id,
        "lithology": part["lithology"]
    }

    for col in ["Pm", "tm", "Tm", "time_interval", "ti", "tf"]:
        first_valid_idx = lowgrad[col].first_valid_index()
        if first_valid_idx is not None:  # Check if a valid index exists
            particle_data[col] = lowgrad.at[first_valid_idx, col]
        else:
            particle_data[col] = np.nan  # No valid value in the column



    if id % 10 == 0:
        fig = plt.figure(figsize=(8, 6))
        # plt.plot(lowgrad["time"], lowgrad["Pm"], color="green", linewidth=2, zorder=5)
        # Plot the time interval calculated before, from ti to tf
        for i, row in lowgrad.iterrows():
            plt.plot([row["ti"], row["tf"]], [row["Pm"], row["Pm"]], color="green", linewidth=2, zorder=5)
        plt.plot(data["time"], data["Plith"], color="grey")
        plt.scatter(lowgrad["tm"], lowgrad["Pm"], color="red", zorder=10)
        plt.xlabel("Time (Myr)")
        plt.ylabel("Pressure (GPa)")
        plt.title(f"Particle {id} - Lithology: {part['lithology']}")
        plt.ylim(0, 3.)
        plt.savefig(f"{sloc}/stagnant_times_{id}.png")
        plt.close()

    lowgrad.to_csv(f"{sfiles}/exhumed_times_{id}.txt", sep=" ", index=False, float_format='%.2f', na_rep="NaN")

    return particle_data


############### MAIN ####################

def main():
    parser = argparse.ArgumentParser(description='Script that gets some models and gives the kinematic indicators')
    parser.add_argument('json_file', help='json file with model name, time at the end of computation, folder where to save the plots, legend')
    args = parser.parse_args()

    json_loc = '/home/vturino/PhD/projects/exhumation/pyInput/'

    with open(f"{json_loc}{args.json_file}") as json_file:
        configs = json.load(json_file)
    m = configs["models"][0]

    plot_loc = f"/home/vturino/PhD/projects/exhumation/plots/single_models/{m}"
    txt_loc = f'{plot_loc}/txt_files'
    if not os.path.exists(txt_loc):
        os.mkdir(txt_loc)

    sloc = f"{plot_loc}/exhumed_times"
    if not os.path.exists(sloc):
        os.mkdir(sloc)

    sfiles = f"{plot_loc}/txt_files/exhumed_times"
    if not os.path.exists(sfiles):
        os.mkdir(sfiles)

    part = pd.read_csv(f"{txt_loc}/exhumed_particles.txt", sep="\s+")
    columns = ["id", "lithology", "Pm", "tm", "Tm", "time_interval", "ti", "tf"]
    particles = pd.DataFrame(columns=columns, index=part.index)

    stagnation_min = 10.

    grad_thresh = 0.008
    time_thresh = 0.5

    with ProcessPoolExecutor(max_workers=2) as executor:
        futures = [
            executor.submit(process_particle, part["id"].iloc[i], part.iloc[i], txt_loc, stagnation_min, time_thresh, sloc, sfiles, grad_thresh)
            for i in range(len(part))
        ]

        for future in tqdm(as_completed(futures), total=len(futures)):
            result = future.result()
            if result:
                for key, value in result.items():
                    particles.at[result["id"], key] = value
    particles.dropna(subset=["Pm", "tm", "Tm", "time_interval", "ti", "tf"], how='all', inplace=True)
    particles = particles.astype({col: 'float' for col in columns if col not in ["id", "lithology"]})
    particles.to_csv(f"{txt_loc}/exhumed_times.txt", sep=" ", index=False, float_format='%.2f', na_rep="NaN")

    print("Initial number of particles: ", len(part))
    print("Exhumed particles after processing: ", len(particles))


if __name__ == "__main__":
    main()