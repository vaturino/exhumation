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
    start_time = None
    for i, row in data.iterrows():
        # If duration == time_thresh, we're still in the same interval
        if row['duration'] == time_thresh:
            if start_time is None:
                start_time = row['time']  # Mark the start time of the interval
            end_time = row['time']  # Update the end time (current row's time)

            # If we reach the last row of consecutive 0.5 intervals or the last row in the bin, assign the time_bin
            if i == len(data) - 1 or data.at[i+1, 'duration'] != time_thresh:
                # Assign the time bin to all rows in the current interval
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
    
def process_times_for_particle(data, stagnation_min, time_thresh):
    data['gradient'] = np.gradient(data["Plith"], data["time"])
    lowgrad = data[abs(data["gradient"]) < 0.01].reset_index(drop=True)
    lowgrad["duration"] = lowgrad["time"].diff()
    lowgrad['time_bin'] = None

    lowgrad = compute_time_intervals(lowgrad, stagnation_min, time_thresh)
    lowgrad = calculate_middle_values(lowgrad)
    return lowgrad








############### MAIN ####################

def main():
    parser = argparse.ArgumentParser(description='Script that gets some models and gives the kinematic indicators')
    parser.add_argument('json_file', help='json file with model name, time at the end of computation, folder where to save the plots, legend')
    args = parser.parse_args()

    json_loc = '/home/vturino/PhD/projects/exhumation/pyInput/'

    # Read the json file
    with open(f"{json_loc}{args.json_file}") as json_file:
        configs = json.load(json_file)
    m = configs["models"][0]

    # Create the folders to save the plots
    plot_loc = f"/home/vturino/PhD/projects/exhumation/plots/single_models/{m}"
    txt_loc = f'{plot_loc}/txt_files'
    if not os.path.exists(txt_loc):
        os.mkdir(txt_loc)

    sloc = f"{plot_loc}/stagnant_times"
    if not os.path.exists(sloc):
        os.mkdir(sloc)

    sfiles = f"{plot_loc}/txt_files/stagnant_times"
    if not os.path.exists(sfiles):
        os.mkdir(sfiles)

    timing = ["dyn", "trans", "kin"]



    fixed_columns = ["id", "lithology"]
    columns = fixed_columns + [f"{c}_{t}" for c in ["Pm", "tm", "Tm", "time_interval", "ti", "tf"] for t in timing]
    
    part = pd.read_csv(f"{txt_loc}/stagnant_particles.txt", sep="\s+")
    particles = pd.DataFrame(columns=columns, index=part.index)
    

    stagnation_min = 10.
    grad_thresh = 0.01
    time_thresh = 0.5

    count = 0

    for i in range(len(part)):
        id = part["id"].iloc[i]
        data = load_data(id, txt_loc)
        data['gradient'] = np.gradient(data["Plith"], data["time"])

        lowgrad = process_times_for_particle(data, stagnation_min, time_thresh)
        lowgrad["lithology"] = part["lithology"].iloc[i]

        bin_number = lowgrad["time_bin"].nunique()
        if bin_number == 0:
            continue

        count += 1


        c = ["Pm", "tm", "Tm", "time_interval", "ti", "tf"]

        # Pre-filter the lowgrad DataFrame once for each condition
        lowgrad_dyn = lowgrad[lowgrad["tm"] < 33.]
        lowgrad_kin = lowgrad[lowgrad["tm"] > 37.]
        lowgrad_trans = lowgrad[(lowgrad["tm"] > 33.) & (lowgrad["tm"] < 37.)]

        # Iterate over each column in `c` and check the subsets
        for col in c:
            # Check and assign for dyn
            if not lowgrad_dyn[col].empty:
                particles[f"{col}_dyn"].iloc[i] = lowgrad_dyn[col].values[0]
            else:
                particles[f"{col}_dyn"].iloc[i] = np.nan

            # Check and assign for kin
            if not lowgrad_kin[col].empty:
                particles[f"{col}_kin"].iloc[i] = lowgrad_kin[col].values[0]
            else:
                particles[f"{col}_kin"].iloc[i] = np.nan

            # Check and assign for trans
            if not lowgrad_trans[col].empty:
                particles[f"{col}_trans"].iloc[i] = lowgrad_trans[col].values[0]
            else:
                particles[f"{col}_trans"].iloc[i] = np.nan

        # if all values of particles are nan, then remove the particle
        if particles.iloc[i].isnull().all():
            part.drop(i, inplace=True)
            continue


        particles["id"].iloc[i] = id
        particles["lithology"].iloc[i] = part["lithology"].iloc[i]

        
        if i % 10 == 0:
            fig = plt.figure(figsize=(8, 6))
            plt.plot(lowgrad["time"], lowgrad["Pm"], color="green", linewidth = 2, zorder = 5)
            plt.plot(data["time"], data["Plith"], color="grey")
            plt.scatter(lowgrad["tm"], lowgrad["Pm"], color="red", zorder = 10)
            plt.xlabel("Time (Myr)")
            plt.ylabel("Pressure (GPa)")
            plt.title(f"Particle {id} - Lithology: {part['lithology'].iloc[i]}")
            plt.ylim(0,2.5)
            plt.savefig(f"{sloc}/stagnant_times_{id}.png")
            plt.close()

        lowgrad.to_csv(f"{sfiles}/stagnant_times_{id}.txt", sep=" ", index=False, float_format='%.2f')

    print(particles)
    print("Number of particles: ", len(part))
    print("Number of particles that stagnate at least once: ", count)   





if __name__ == "__main__":
    main()

