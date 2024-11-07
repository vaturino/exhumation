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
import math as math
from scipy.signal import savgol_filter 
import seaborn as sns
from matplotlib.axes._axes import _log as matplotlib_axes_logger
matplotlib_axes_logger.setLevel('ERROR')
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor

path = str(Path(Path(__file__).parent.absolute()).parent.absolute())
sys.path.insert(0, path)
from libraries.functions import *
from libraries.particles import *

def parse_arguments():
    parser = argparse.ArgumentParser(description='Script that gets some models and gives the kinematic indicators')
    parser.add_argument('json_file', help='json file with model name, time at the end of computation, folder where to save the plots, legend')
    return parser.parse_args()

def load_configs(json_loc, json_file):
    with open(f"{json_loc}{json_file}") as file:
        return json.load(file)

def setup_directories(plot_loc):
    txt_loc = f'{plot_loc}/txt_files'
    pt_files = f'{txt_loc}/PT'
    os.makedirs(pt_files, exist_ok=True)
    return txt_loc, pt_files

def initialize_particle_files(pt_files, npart, compositions):
    for i in range(npart):
        with open(f"{pt_files}/pt_part_{i}.txt", "w+") as pt:
            pt.write("id time x y P Plith T depth vx vy " + " ".join(compositions) + "\n")

def process_time_step(t, csvs_loc, m, compositions, tr, idx, cmyr, rhog, ymax=901.e3):
    fcol = ["id", "position:0", "position:1", "p", "T", "position:1", "velocity:0", "velocity:1"]
    df = pd.read_parquet(f"{csvs_loc}{m}/particles/full.{t}.gzip", columns=fcol + compositions)
    fil = (df[compositions] > tr).any(axis=1)
    df = df[fil]
    conds = df[df['id'].isin(idx['id'])].sort_values(by='id')

    data = []
    for i in range(len(idx)):
        row = [
            conds["id"].iloc[i], t, conds["position:0"].iloc[i], conds["position:1"].iloc[i],
            conds["p"].iloc[i] / 1e9, (ymax - conds["position:1"].iloc[i])*rhog/1.e9, (conds["T"].iloc[i] + 0.6*(ymax - conds["position:1"].iloc[i])/1.e3) - 273.,
            conds["position:1"].iloc[i] / 1.e3,
            conds["velocity:0"].iloc[i] * cmyr,
            conds["velocity:1"].iloc[i] * cmyr
        ]
        row.extend(conds[col].iloc[i] for col in compositions)
        data.append(row)
    return data

def main():
    args = parse_arguments()
    json_loc = '/home/vturino/PhD/projects/exhumation/pyInput/'
    models_loc = '/home/vturino/Vale_nas/exhumation/raw_outputs/'
    csvs_loc = '/home/vturino/Vale_nas/exhumation/gz_outputs/'
    cmyr = 1e2 * (60 * 60 * 24 * 365)
    tr = 1.e-2
    rho = 3100
    g = 9.81
    rhog = rho * g
    ymax = 901.e3

    configs = load_configs(json_loc, args.json_file)
    compositions = configs['compositions']

    for ind_m, m in tqdm(enumerate(configs['models'])):
        time_array = np.zeros((len(os.listdir(f"{csvs_loc}{m}/particles")), 2))
        stat = pd.read_csv(f"{models_loc}{m}/statistics", skiprows=configs['head_lines'], sep='\s+', header=None)
        time_array = grab_dimTime_particles(f"{csvs_loc}{m}/particles", stat, time_array, configs['head_lines'] - 1)

        plot_loc = f"/home/vturino/PhD/projects/exhumation/plots/single_models/{configs['models'][ind_m]}"
        txt_loc, pt_files = setup_directories(plot_loc)

        with open(f"{txt_loc}/particles_indexes.txt", "r") as f:
            npart = len(f.readlines()) - 1

        initialize_particle_files(pt_files, npart, compositions)
        idx = pd.read_csv(f"{txt_loc}/particles_indexes.txt", sep="\s+")

        all_data = [[] for _ in range(npart)]

        with ProcessPoolExecutor() as executor:
            futures = [executor.submit(process_time_step, t, csvs_loc, m, compositions, tr, idx, cmyr, rhog, ymax) for t in range(1, len(time_array))]
            for future in tqdm(futures):
                data = future.result()
                for i, row in enumerate(data):
                    all_data[i].append(row)

        for i in range(npart):
            all_data[i].sort(key=lambda x: x[1])  # Sort by time
            with open(f"{pt_files}/pt_part_{i}.txt", "a+") as pt:
                for row in all_data[i]:
                    pt.write(" ".join(f"{val:.3f}" if isinstance(val, float) else f"{val:.0f}" for val in row) + "\n")

if __name__ == "__main__":
    main()
