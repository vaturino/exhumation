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


def load_data(id, txt_loc):
    data = pd.read_csv(f"{txt_loc}/PT/pt_part_{id}.txt", sep="\s+")
    data["time"] = data["time"] / 2
    data["ts"] = data["time"] * 2
    return data

def avg_pressures_and_extremes(data, time_interval):
    tf = data["time"].max()
    ts = tf - time_interval
    Pf = data[data["time"] == tf]["Plith"].values[0]
    Ps = data[data["time"] == ts]["Plith"].values[0]
    Pm = (Pf + Ps) / 2
    return Pm, Pf, Ps, tf, ts

def compute_time_intervals(data, Pthresh, Pm):
    index = []
    for k in range(len(data["Plith"]) - 1, -1, -1):
        if abs(data["Plith"].iloc[k] - Pm) <= Pthresh:
            index.append(k)
    
    for j in range(len(index) - 1):
        if index[j] - index[j + 1] != 1:
            index = index[:j + 1]
            break

    time_interval = pd.DataFrame(data["time"].iloc[index[-1]:index[0]])
    time_interval["Pm"] = Pm

    return time_interval


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

    
    part = pd.read_csv(f"{txt_loc}/stagnant_particles.txt", sep="\s+")
    particles = pd.DataFrame(columns=["id", "lithology", "Pm", "tm", "time_interval"], index=range(len(part)))

    for i in range(len(part)):
        id = part["id"].iloc[i]
        data = load_data(id, txt_loc)

        thresh = min(data["Plith"].max() / 20, 0.05)
        time_interval = 15
        Pm, Pf, Ps, tf, ts = avg_pressures_and_extremes(data, time_interval)

        time_interval = compute_time_intervals(data, thresh, Pm)

        if time_interval["time"].max() - time_interval["time"].min() > 10:
            particles.at[i, "id"] = id
            particles.at[i, "lithology"] = part["lithology"].iloc[i]
            particles.at[i, "Pm"] = Pm
            particles.at[i, "tm"] = (tf + ts) / 2
            particles.at[i, "time_interval"] = time_interval["time"].max() - time_interval["time"].min()
        else:
            continue

        if i % 100 == 0:
            plt.plot(time_interval["time"], time_interval["Pm"], color="green")
            plt.plot(data["time"], data["Plith"], color="grey")
            plt.scatter(tf, Pf, color="blue")
            plt.scatter(ts, Ps, color="blue")
            plt.scatter((tf + ts) / 2, Pm, color="red")
            plt.axhline(y=Pm + thresh, color='r', linestyle='--')
            plt.axhline(y=Pm - thresh, color='r', linestyle='--')
            plt.xlabel("Time (Myr)")
            plt.ylabel("Pressure (GPa)")
            plt.title(f"Particle {id} - Lithology: {part['lithology'].iloc[i]}")
            plt.ylim(0, data["Plith"].max() + thresh)
            plt.savefig(f"{sloc}/stagnant_times_{id}.png")
            plt.close()

    particles.dropna(subset=["id"], inplace=True)
    float_cols = particles.columns.difference(["id", "lithology"])
    particles[float_cols] = particles[float_cols].astype(float)
    particles.to_csv(f"{txt_loc}/stagnant_times.txt", sep="\t", index=False, header=True, float_format='%.2f')


if __name__ == "__main__":
    main()



# def compute_time_intervals(data, Pthresh, min_delta_time=10):
#     intervals = []
#     start, end = None, None
    
#     for i in range(len(data) - 1):
#         if data["time"].iloc[i] > 10 and abs(data["Plith"].iloc[i] - data["Plith"].iloc[i + 1]) <= Pthresh:
#             if start is None:
#                 start = i
#             end = i + 1
#         else:
#             if start is not None and (data["time"].iloc[end] - data["time"].iloc[start]) >= min_delta_time:
#                 Pm = data["Plith"].iloc[start:end + 1].mean()
#                 interval = {
#                     "start_time": data["time"].iloc[start],
#                     "end_time": data["time"].iloc[end],
#                     "Pm": Pm,
#                     "duration": data["time"].iloc[end] - data["time"].iloc[start]
#                 }
#                 intervals.append(interval)
#             start, end = None, None

#     if start is not None and (data["time"].iloc[end] - data["time"].iloc[start]) >= min_delta_time:
#         Pm = data["Plith"].iloc[start:end + 1].mean()
#         interval = {
#             "start_time": data["time"].iloc[start],
#             "end_time": data["time"].iloc[end],
#             "Pm": Pm,
#             "duration": data["time"].iloc[end] - data["time"].iloc[start]
#         }
#         intervals.append(interval)

#     return intervals


        # intervals = compute_time_intervals(data, thresh)
        # rows_to_add = []

        # for interval in intervals:
        #     if interval["duration"] >= 10:
        #         rows_to_add.append({
        #             "id": id,
        #             "lithology": part["lithology"].iloc[i],
        #             "Pm": interval["Pm"],
        #             "tm": (interval["start_time"] + interval["end_time"]) / 2,
        #             "time_interval": interval["duration"]
        #         })

        # if rows_to_add:
        #     particles = pd.concat([particles, pd.DataFrame(rows_to_add)], ignore_index=True)