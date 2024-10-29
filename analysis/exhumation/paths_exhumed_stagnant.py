#! /usr/bin/python3

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

def plot_particle(case, p, elist, slist, pt_files, ymax, eloc, name, location):
    id = eval(case)['id'].iloc[p]
    pt_single = pd.read_csv(f"{pt_files}/pt_part_{id}.txt", sep="\s+")
    pt_single["Plith"] = (ymax - pt_single["depth"]) * 1e3 * 9.81 * 3100 / 1e9
    eval(case)["tin"] = 0

    k = -100

    for i in range(len(pt_single)):
        if pt_single.depth.iloc[0] - pt_single.depth.iat[i] >= 2.:
            eval(case)["tin"].iloc[p] = pt_single["time"].iat[i]/2.
            k = i
            break

    fig, ax = plt.subplots(1, 2, figsize=(12, 6))
    ax[0].plot(pt_single["T"] + 0.6 * (pt_single["depth"].max() - pt_single["depth"]) - 273, pt_single["Plith"], label='Pressure')
    # ax[0].plot(pt_single["T"] + 0.6 * (pt_single["depth"].max() - pt_single["depth"]) - 273, pt_single["Plith"].rolling(window=2).mean(), label='Smoothed Pressure', linestyle='--')
    ax[0].scatter(eval(case)['maxPT'].iloc[p], eval(case)['maxPP'].iloc[p], color='red')
    ax[0].annotate(f"tmax = {eval(case)['tmax'].iloc[p]} Myr", (eval(case)['maxPT'].iloc[p], eval(case)['maxPP'].iloc[p]))
    ax[0].scatter(pt_single["T"].iloc[k]+ 0.6*(pt_single["depth"].max() - pt_single["depth"].iloc[k])-273, pt_single["Plith"].iloc[k], color='green')
    ax[0].annotate(f"tin = {eval(case)['tin'].iloc[p]} Myr", (pt_single["T"].iloc[k]-273, pt_single["Plith"].iloc[k]))
    ax[0].set_xlabel("Temperature (°C)")
    ax[0].set_ylabel("Pressure (GPa)")

    ax[1].plot(pt_single["time"], pt_single["Plith"], label='Trajectory')
    ax[1].set_xlabel("Time (Myr)")
    ax[1].set_ylabel("Pressure (GPa)")
    ax[1].invert_yaxis()
    plt.suptitle(f"Particle {id} - {name} lithology: {eval(case)['lithology'].iloc[p]}")
    plt.savefig(f"{location}/{name}_particle_{id}.png", dpi=500)
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

    sloc = f"{plot_loc}/stagnant"
    os.makedirs(sloc, exist_ok=True)

    pt_files = f'{txt_loc}/PT'

    elist = pd.read_csv(f"{txt_loc}/exhumed_particles.txt", sep="\s+")
    slist = pd.read_csv(f"{txt_loc}/stagnant_particles.txt", sep="\s+")

    # max plot (no free-surface)
    ymax = 900.

    cases = ["elist", "slist"]
    names = ["Exhumed", "Stagnant"]
    locations = [eloc, sloc]

    with ProcessPoolExecutor() as executor:
        futures = []
        for case, name, location in zip(cases, names, locations):
            for p in range(0, len(eval(case)), 100):
                futures.append(executor.submit(plot_particle, case, p, elist, slist, pt_files, ymax, eloc, name, location))

        for future in tqdm(as_completed(futures), total=len(futures)):
            future.result()

if __name__ == "__main__":
    main()
