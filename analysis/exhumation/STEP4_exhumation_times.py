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
from scipy.stats import linregress


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
    texh_thresh = 0.35

    n_old = len(part)
    n_new = 0

    exhumed_times = pd.DataFrame(index=range(len(part)), columns=["id", "lithology", "tin", "Pin", "tmax", "Pmax", "texh", "Pexh", "tfin", "Pfin"])

    for i in tqdm(range(0, len(part))):
        id = part["id"].iloc[i]
        particle = load_data(id, txt_loc)
        # particle["Plith"] = particle["Plith"].rolling(window=10, min_periods=1).mean()
        # particle["time"] = particle["time"].rolling(window=10, min_periods=1).mean()

        if particle.empty:
            continue

        Pmax = particle["Plith"].max()
        idmax = particle["Plith"].idxmax()
        Pexh = (1-texh_thresh) * Pmax

        time = 50.
        idxtime = particle[particle["time"] == time].index[0]
        # idxtime = particle.iloc[(particle["time"] - time).abs().argsort()[:1]].index[0]


        if particle["Plith"].iloc[idxtime] <= Pexh:

            partial = particle.loc[idmax:].reset_index(drop=True)
            time = partial["time"].values
            Plith = partial["Plith"].values


            gradient = np.gradient(Plith, time)
            min_duration = 5.0  # Myr

            tin = None
            Pin = None

            for start_idx in range(len(time) - 2):
                t_start = time[start_idx]
                idx = start_idx

                while idx + 1 < len(time) and gradient[idx] < 0:
                    idx += 1
                    if time[idx] - t_start >= min_duration:
                        tin = t_start + 1 # Add 1 Myr to avoid the exact Pmax time
                        Pin = Plith[np.isclose(time, tin)][0]
                        n_new += 1
                        break

                if tin is not None:
                    break
            
            
            texh = partial[partial["Plith"] <= Pexh]["time"].min()
            tfin = particle["time"].iloc[-1]
            Pfin = particle["Plith"].iloc[-1]
            lithology = part["lithology"].iloc[i]

            # Plotting
            if i % 50 == 0:
                plt.plot(particle["time"], particle["Plith"], alpha=0.5)
                plt.axhline(y=Pexh, color='red', linestyle='--', label='Exhumation Threshold')
                plt.axvline(x=50, color='grey', linestyle='-', linewidth=1, alpha = 0.5)
                if tin is not None:
                    plt.scatter(tin, Pin, color='blue', label='Tin', s=50)
                plt.legend()
                plt.xlabel("Time (Myr)")
                plt.ylabel("Pressure (GPa)")
                plt.title(f"Particle {id}: {part['lithology'].iloc[i]}")
                plt.savefig(f"{sloc}/particle_{id}.png")
                plt.clf()


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
    

    print(f"Old particles: {n_old}, New particles: {n_new}")



if __name__ == "__main__":
    main()
