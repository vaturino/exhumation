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


import pandas as pd
import numpy as np

def load_data(id, txt_loc):
    data = pd.read_csv(f"{txt_loc}/PT/pt_part_{id}.txt", sep="\s+")
    data["time"] = data["time"] / 2
    data["ts"] = data["time"] * 2
    return data

def compute_stagnant_maxima(data, Pthresh):
    # Calculate Pdiff using groupby once
    data["Pdiff"] = data.groupby("time_bin")["Plith"].transform(lambda x: abs(x.max() - x.min()))
    data.loc[data["time"] <= 10., "Pdiff"] = np.nan
    Pmaxloc_medians = data[data["Pdiff"] <= Pthresh].groupby("time_bin")["Plith"].median()
    data["Pmaxloc"] = data["time_bin"].map(Pmaxloc_medians)

    # Step 1: Mark non-NaN Pmaxloc values
    data['Pmaxloc_not_nan'] = data['Pmaxloc'].notna()

    # Step 2: Calculate time_bin_consecutive without loops
    data['time_bin_consecutive'] = np.where(
        data['Pmaxloc_not_nan'],
        (data['Pmaxloc_not_nan'] & (data['time_bin'] != data['time_bin'].shift())).cumsum(),
        np.nan
    )
    
    # Conditionally update time_bin_consecutive with minimal looping
    mask = data['time_bin_consecutive'].notna() & (data['time_bin_consecutive'] > data['time_bin_consecutive'].shift())
    data.loc[mask, 'time_bin_consecutive'] = np.where(
        abs(data.loc[mask, "Pmaxloc"].values - data.loc[mask, "Pmaxloc"].shift().values) < Pthresh,
        data.loc[mask, 'time_bin_consecutive'].shift().values,
        data.loc[mask, 'time_bin_consecutive'].values
    )

    # Determine particle direction based on last 10 Plith differences
    lastPdata = data["Plith"].iloc[-1] - data["Plith"].iloc[-10]
    data["direction"] = np.select(
        [lastPdata <= -data["Plith"].max() / 100., lastPdata > data["Plith"].max() / 100.],
        ["rising", "subducting"],
        default="stagnant"
    )
    
    # Calculate Pdiff_new and Pmaxloc_new using optimized groupby
    data['Pdiff_new'] = data.groupby('time_bin_consecutive')["Plith"].transform(lambda x: x.iloc[-1] - x.iloc[0])
    data['Pmaxloc_new'] = data.groupby('time_bin_consecutive')['Pmaxloc'].transform("mean")

    return data

def process_particle(p, part, txt_loc):
    id = part["id"].iloc[p]
    data = load_data(id, txt_loc)

    Pthresh = data["Plith"].max() / 20

    # Bin data by time intervals of 10 Myr
    bins = np.arange(0, 51, 10)
    data["time_bin"] = pd.cut(data["time"], bins)

    data = compute_stagnant_maxima(data, Pthresh)

    # Check if the last valid Pmaxloc_new value is in the last time_bin
    last_valid_index = data["Pmaxloc_new"].last_valid_index()
    if last_valid_index is None or data["time_bin"].iloc[last_valid_index] != data["time_bin"].max():
        return None, None  # Particle does not meet the condition

    # Collect particle data
    particle_data = {
        "id": id,
        "lithology": part["lithology"].iloc[p],
        "direction": data["direction"].iloc[0],
        "timePmax": data["time"].iloc[data["Plith"].idxmax()],
        "Pmax": data["Plith"].max(),
        "Tmax": data["T"].iloc[data["Plith"].idxmax()],
    }

    # Populate data for unique Pmaxloc_new values
    unique_pmaxloc_new = data["Pmaxloc_new"].dropna().unique()
    unique_consecutives = data['time_bin_consecutive'].dropna().unique()
    for k, (pmax, consec) in enumerate(zip(unique_pmaxloc_new, unique_consecutives), start=1):
        particle_data[f"P{k}"] = pmax
        particle_data[f"timeP{k}"] = data.loc[data["time_bin_consecutive"] == consec, "time"].median()
        particle_data[f"T{k}"] = data.loc[data["time_bin_consecutive"] == consec, "T"].mean()

    return particle_data, data




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

    sloc = f"{plot_loc}/stagnant_maxima"
    if not os.path.exists(sloc):
        os.mkdir(sloc)

    part = pd.read_csv(f"{txt_loc}/stagnant_particles.txt", sep="\s+")
    particles = pd.DataFrame(columns=["id", "lithology",  "direction", "timePmax", "Pmax", "Tmax", "timeP1", "P1", "T1", "timeP2", "P2", "T2", "timeP3", "P3", "T3", "timeP4", "P4", "T4",], index=range(len(part)))

    with ProcessPoolExecutor() as executor:
        futures = {executor.submit(process_particle, p, part, txt_loc): p for p in range(len(part))}
        # Inside the main function, after obtaining particle_data and data
        for future in tqdm(as_completed(futures), total=len(futures)):
            p = futures[future]
            particle_data, data = future.result()

            # Skip particles that do not meet the condition
            if particle_data is None:
                continue

            for key, value in particle_data.items():
                particles.at[p, key] = value

            if p % 100 == 0:
                plt.plot(data["time"], data["Plith"], color="grey")
                plt.plot(data["time"], data["Pmaxloc_new"], color="blue")
                plt.xlabel("Time (Myr)")
                plt.ylabel("Pressure (GPa)")
                plt.title(f"Particle {particle_data['id']} - Lithology: {particle_data['lithology']} - Direction: {particle_data['direction']}")
                plt.ylim(0, data["Plith"].max() + 0.1)
                plt.savefig(f"{sloc}/stagnant_maxima_{particle_data['id']}.png")
                plt.close()


            # Save only the filtered particles
            particles.dropna(subset=["id"], inplace=True)
            particles.to_csv(f"{txt_loc}/stagnant_maxima.txt", sep="\t", index=False, header=True, float_format='%.2f')

    print(f"Number of actually stagnant particles: {len(particles)}")

if __name__ == "__main__":
    main()















# def process_particle(p, part, txt_loc):
#     id = part["id"].iloc[p]
#     data = load_data(id, txt_loc)

#     Pthresh = data["Plith"].max()/20

#     # Bin data by time intervals of 10 Myr
#     bins = np.arange(0, 51, 10)
#     data["time_bin"] = pd.cut(data["time"], bins)

#     data = compute_stagnant_maxima(data, Pthresh)

#     particle_data = {
#         "id": id,
#         "lithology": part["lithology"].iloc[p],
#         "direction": data["direction"].iloc[0],
#         "timePmax": data["time"][data["Plith"].idxmax()],
#         "Pmax": data["Plith"].max(),
#         "Tmax": data["T"][data["Plith"].idxmax()],
#     }

#     for k in range(len(data["Pmaxloc_new"].dropna().unique())): 
#         particle_data[f"P{k+1}"] = data["Pmaxloc_new"].dropna().unique()[k]
#         particle_data[f"timeP{k+1}"] = data[data["time_bin_consecutive"] == data["time_bin_consecutive"].dropna().unique()[k]]["time"].median()
#         particle_data[f"T{k+1}"] = data[data["time_bin_consecutive"] == data["time_bin_consecutive"].dropna().unique()[k]]["T"].mean()

#     return particle_data, data

# def main():
#     parser = argparse.ArgumentParser(description= 'Script that gets some models and gives the kinematic indicators')
#     parser.add_argument('json_file', help='json file with model name, time at the end of computation, folder where to save the plots, legend')
#     args = parser.parse_args()

#     json_loc = '/home/vturino/PhD/projects/exhumation/pyInput/'

#     # Read the json file
#     with open(f"{json_loc}{args.json_file}") as json_file:
#             configs = json.load(json_file)
#     m = configs["models"][0]

#     # Create the folders to save the plots
#     plot_loc = f"/home/vturino/PhD/projects/exhumation/plots/single_models/{m}"
#     txt_loc = f'{plot_loc}/txt_files'
#     if not os.path.exists(txt_loc):
#         os.mkdir(txt_loc)

#     sloc = f"{plot_loc}/stagnant_maxima"
#     if not os.path.exists(sloc):
#         os.mkdir(sloc)

#     part = pd.read_csv(f"{txt_loc}/stagnant_particles.txt", sep="\s+")
#     particles = pd.DataFrame(columns=["id", "lithology",  "direction", "timePmax", "Pmax", "Tmax", "timeP1", "P1", "T1", "timeP2", "P2", "T2", "timeP3", "P3", "T3", "timeP4", "P4", "T4",], index=range(len(part)))

#     with ProcessPoolExecutor() as executor:
#         futures = {executor.submit(process_particle, p, part, txt_loc): p for p in range(len(part))}
#         for future in tqdm(as_completed(futures), total=len(futures)):
#             p = futures[future]
#             particle_data, data = future.result()

#             for key, value in particle_data.items():
#                 particles.at[p, key] = value

#             if p % 10 == 0:
#                 plt.plot(data["time"], data["Plith"], color="grey")
#                 plt.plot(data["time"], data["Pmaxloc_new"], color="blue")
#                 plt.xlabel("Time (Myr)")
#                 plt.ylabel("Pressure (GPa)")
#                 plt.title(f"Particle {particle_data['id']} - Lithology: {particle_data['lithology']} - Direction: {particle_data['direction']}")
#                 plt.ylim(0, 2.0)
#                 plt.savefig(f"{sloc}/stagnant_maxima_{particle_data['id']}.png")
#                 plt.close()


#     particles.to_csv(f"{txt_loc}/stagnant_maxima.txt", sep="\t", index=False, header=True, float_format='%.2f')

# if __name__ == "__main__":
#     main()