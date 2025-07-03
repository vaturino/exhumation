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
    ti_thresh = 0.20
    texh_thresh = 0.35

    exhumed_times = pd.DataFrame(index=range(len(part)), columns=["id", "lithology", "tin", "Pin", "tmax", "Pmax", "texh", "Pexh", "tfin", "Pfin"])

    for i in tqdm(range(0, len(part))):
        id = part["id"].iloc[i]
        particle = load_data(id, txt_loc)

        if particle.empty:
            continue
            
        Pmax = particle["Plith"].max()
        idmax = particle["Plith"].idxmax()
        partial = particle.loc[idmax:]

        Pin = (1-ti_thresh)* Pmax
        tin_candidates = partial[partial["Plith"] <= Pin]
        tin = None
        for idx in tin_candidates.index:
            if idx + 5 < len(particle) and all(particle["Plith"].iloc[idx:idx+5].diff().dropna() < 0):
                tin = tin_candidates.loc[idx, "time"]
            break
        Pexh = (1-texh_thresh) * Pmax
        texh = partial[partial["Plith"] <= Pexh]["time"].min()
        tfin = particle["time"].iloc[-1]
        Pfin = particle["Plith"].iloc[-1]
        lithology = part["lithology"].iloc[i]

        exhumed_times.iloc[i] = {
            "id": id,
            "tin": tin,
            "Pin": round(Pin, 3),
            "tmax": particle["time"].iloc[idmax],
            "Pmax": round(Pmax, 3),
            "texh": texh,
            "Pexh": round(Pexh, 3),
            "tfin": tfin,
            "Pfin": round(Pfin, 3),
            "lithology": lithology
        }
    
    exhumed_times.to_csv(f"{txt_loc}/exhumed_times.txt", sep="\t", index=False, float_format='%.2f', na_rep="NaN")


if __name__ == "__main__":
    main()
