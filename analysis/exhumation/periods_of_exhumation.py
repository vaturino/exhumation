#! /usr/bin/python3

import pandas as pd
import numpy as np
import json
from tqdm import tqdm
import matplotlib.pyplot as plt
import argparse
import os
import warnings
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

warnings.simplefilter(action='ignore', category=FutureWarning)
pd.options.mode.chained_assignment = None  # default='warn'

path = str(Path(Path(__file__).parent.absolute()).parent.absolute())
import sys
sys.path.insert(0, path)
from libraries.functions import *
from libraries.particles import *
from libraries.exhumation import *  

############### FUNCTIONS ####################

def load_data(id, txt_loc):
    data = pd.read_csv(f"{txt_loc}/PT/pt_part_{id}.txt", sep="\s+")
    data["time"] = data["time"] / 2
    data["ts"] = data["time"] * 2
    return data

def compute_time_intervals(data, stagnation_min, time_thresh, Pthresh):
    texh_candidates = data[(data["Plith"] <= Pthresh) & (data["time"] >= stagnation_min)]
    texh = texh_candidates["time"].min() if not texh_candidates.empty else data["time"].max() - 0.5

    data['gradient'] = np.gradient(data["Plith"], data["time"])
    start_time = None

    for i, row in data.iterrows():
        if row['duration'] == time_thresh:
            if start_time is None:
                start_time = row['time']
            end_time = row['time']
            window = data.loc[i:min(i + 4, len(data)-1), 'gradient']
            decreasing_grad = (window < 0.1).mean() >= 0.8
            if (i == len(data) - 1 or data.at[i + 1, 'duration'] != time_thresh) and decreasing_grad:
                for j in range(i, -1, -1):
                    if data.at[j, 'duration'] == time_thresh:
                        adjusted_start_time = max(start_time, texh)
                        if adjusted_start_time < end_time:
                            data.at[j, 'time_bin'] = f"[{adjusted_start_time}, {end_time})"
                            data.at[j, 'ti'] = adjusted_start_time
                            data.at[j, 'tf'] = end_time
                            data.at[j, 'time_interval'] = end_time - adjusted_start_time
                    else:
                        break
                start_time = None
        else:
            start_time = None
    return data

def calculate_middle_values(data):
    avg_values = data.groupby('time_bin').agg({'Plith': 'mean', 'time': 'mean', 'T': 'mean'}).reset_index()
    data["Pm"] = np.nan
    data["tm"] = np.nan
    data["Tm"] = np.nan
    for tint in avg_values["time_bin"]:
        mask = data["time_bin"] == tint
        data.loc[mask, "Pm"] = avg_values.loc[avg_values["time_bin"] == tint, "Plith"].values[0]
        data.loc[mask, "tm"] = avg_values.loc[avg_values["time_bin"] == tint, "time"].values[0]
        data.loc[mask, "Tm"] = avg_values.loc[avg_values["time_bin"] == tint, "T"].values[0]
    return data

def process_times_for_particle(data, stagnation_min, time_thresh, grad_thresh, Pfrac, min_time=15.):
    data['gradient'] = np.gradient(data["Plith"], data["time"])
    Pmax = data["Plith"].max()
    Pthresh = round(Pfrac * Pmax, 3)
    fil = ((data["Plith"] <= Pthresh) & (data["time"] > min_time))
    valid_data = data[fil].copy().reset_index(drop=True)
    valid_data = valid_data[(valid_data["gradient"] < 0) | (abs(valid_data["gradient"]) < grad_thresh)].reset_index(drop=True)
    valid_data["duration"] = valid_data["time"].diff()
    valid_data['time_bin'] = None
    valid_data = compute_time_intervals(valid_data, stagnation_min, time_thresh, Pthresh)
    valid_data = calculate_middle_values(valid_data)
    return valid_data, Pthresh

def process_particle(id, part, txt_loc, stagnation_min, time_thresh, sloc, sfiles, grad_thresh, Pfrac):
    data = load_data(id, txt_loc)
    data["Plith"] = data["Plith"].rolling(window=5, min_periods=1).mean()
    data['gradient'] = np.gradient(data["Plith"], data["time"])
    lowgrad, Pthresh = process_times_for_particle(data, stagnation_min, time_thresh, grad_thresh, Pfrac)
    lowgrad["lithology"] = part["lithology"]
    particle_data = {
        "id": id,
        "lithology": part["lithology"],
        "Pm": np.nan,
        "tm": np.nan,
        "Tm": np.nan,
        "time_interval": np.nan,
        "ti": np.nan,
        "tf": np.nan,
        "texh": np.nan,
        "t_thresh": np.nan
    }
    if lowgrad["time_bin"].nunique() > 0:
        for col in ["Pm", "tm", "Tm", "time_interval", "ti", "tf"]:
            val = lowgrad[col].dropna().iloc[0] if not lowgrad[col].dropna().empty else np.nan
            particle_data[col] = val
        particle_data["texh"] = lowgrad["ti"].min()
    t_thresh_candidates = data[data["Plith"] <= Pthresh]
    particle_data["t_thresh"] = t_thresh_candidates["time"].min() if not t_thresh_candidates.empty else np.nan
    if id % 10 == 0:
        fig = plt.figure(figsize=(8, 6))
        for _, row in lowgrad.iterrows():
            plt.plot([row["ti"], row["tf"]], [row["Pm"], row["Pm"]], color="green", linewidth=2)
        plt.plot(data["time"], data["Plith"], color="grey")
        plt.axhline(y=Pthresh, color="lightblue", linestyle="--", linewidth=1)
        plt.scatter(lowgrad["ti"], lowgrad["Pm"], color="red")
        plt.xlabel("Time (Myr)")
        plt.ylabel("Pressure (GPa)")
        plt.title(f"Particle {id} - Lithology: {part['lithology']}")
        plt.ylim(0, 3.)
        plt.savefig(f"{sloc}/stagnant_times_{id}.png")
        plt.close()
    lowgrad.to_csv(f"{sfiles}/exhumed_times_{id}.txt", sep=" ", index=False, float_format='%.2f', na_rep="NaN")
    return particle_data

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
    os.makedirs(txt_loc, exist_ok=True)
    sloc = f"{plot_loc}/exhumed_times"
    os.makedirs(sloc, exist_ok=True)
    sfiles = f"{plot_loc}/txt_files/exhumed_times"
    os.makedirs(sfiles, exist_ok=True)
    part = pd.read_csv(f"{txt_loc}/exhumed_particles.txt", sep="\s+")
    columns = ["id", "lithology", "Pm", "tm", "Tm", "time_interval", "ti", "tf", "texh", "t_thresh"]
    particles = pd.DataFrame(columns=columns, index=part.index)
    stagnation_min = 5.
    grad_thresh = 0.1
    time_thresh = 0.5
    Pfrac = 0.65
    with ProcessPoolExecutor(max_workers=2) as executor:
        futures = [
            executor.submit(process_particle, part["id"].iloc[i], part.iloc[i], txt_loc, stagnation_min, time_thresh, sloc, sfiles, grad_thresh, Pfrac)
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
