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
from scipy.signal import savgol_filter


path = str(Path(Path(__file__).parent.absolute()).parent.absolute())
sys.path.insert(0, path)
from libraries.functions import *
from libraries.particles import *
from libraries.exhumation import *  
import plotly.express as px



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

    sloc = f"{plot_loc}/stagnant_times"
    if not os.path.exists(sloc):
        os.mkdir(sloc)

    ecloc = f"{plot_loc}/stagnant_times_ecl"
    if not os.path.exists(ecloc):
        os.mkdir(ecloc)

    sfiles = f"{plot_loc}/txt_files/stagnant_times"
    if not os.path.exists(sfiles):
        os.mkdir(sfiles)

    timing = ["dyn", "trans", "kin"]

    fixed_columns = ["id", "lithology"]
    columns = fixed_columns + [f"{c}_{t}" for c in ["Pm", "tm", "Tm", "time_interval", "ti", "tf"] for t in timing]
    
    part = pd.read_csv(f"{txt_loc}/stagnant_particles.txt", sep="\s+")

    #find maximum between Pm_dyn, Pm_kin, Pm_trans for all the dataset
    max_Pm = part[["Pm_dyn", "Pm_kin", "Pm_trans"]].max().max() 

    for i in tqdm(range(0, len(part), 10)):
        id = part["id"].iloc[i]
        particle = pd.read_csv(f"{txt_loc}/PT/pt_part_{id}.txt", sep="\s+")
        particle["time"] = particle["time"]/2.
        # particle["Plith"] = particle["Plith"].rolling(window=10, min_periods=1).mean()
        # particle["time"] = particle["time"].rolling(window=10, min_periods=1).mean()

        plt.plot(particle["time"], particle["Plith"], color = "grey", linewidth = 1)
        suff = ["kin", "dyn", "trans"]
        for s in suff:
            plt.scatter(part.iloc[i][f"tm_{s}"], part.iloc[i][f"Pm_{s}"], color="red", zorder = 10)
            plt.hlines(part.iloc[i][f"Pm_{s}"], part.iloc[i][f"ti_{s}"], part.iloc[i][f"tf_{s}"], color="green")

        plt.xlabel("Time (Myr)")
        plt.ylabel("Pressure (GPa)")
        plt.ylim(0, max_Pm+0.1)
        plt.xlim(0, 50)
        plt.title(f"Particle {id} - Lithology: {part.iloc[i]['lithology']}")
        plt.savefig(f"{sloc}/stagnant_times_{id}.png")
        plt.close()

        if part.iloc[i]['lithology'] == "ecl":
            plt.plot(particle["time"], particle["Plith"], color = "grey", linewidth = 1)
            for s in suff:
                plt.scatter(part.iloc[i][f"tm_{s}"], part.iloc[i][f"Pm_{s}"], color="red", zorder = 10)
                plt.hlines(part.iloc[i][f"Pm_{s}"], part.iloc[i][f"ti_{s}"], part.iloc[i][f"tf_{s}"], color="green")

            plt.xlabel("Time (Myr)")
            plt.ylabel("Pressure (GPa)")
            plt.ylim(0, max_Pm+0.1)
            plt.xlim(0, 50)
            plt.title(f"Particle {id} - Lithology: {part.iloc[i]['lithology']}")
            plt.savefig(f"{ecloc}/stagnant_times_{id}.png")
            plt.close()

    # Save the data
    part.to_csv(f"{sfiles}/stagnant_times.txt", sep=" ", columns=columns, index=False)


if __name__ == "__main__":
    main()
