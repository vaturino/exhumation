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
path = str(Path(Path(__file__).parent.absolute()).parent.absolute())
sys.path.insert(0, path)
from libraries.functions import *
from libraries.particles import *
from libraries.exhumation import *  
import plotly.express as px


def load_data(id, txt_loc):
    data = pd.read_csv(f"{txt_loc}/PT/pt_part_{id}.txt", sep="\s+")
    data["Plith"] = (901.e3 - data["y"])*3100*9.81/1e9
    data["time"] = data["time"]/2
    data["ts"] = data["time"]*2
    data["T"] = (data["T"] + 0.6*(data["depth"].max() - data["depth"])) - 273.
    return data

def compute_stagnant_maxima(data, Pthresh):
    data["Pdiff"] = data.groupby("time_bin")["Plith"].transform(lambda x: abs(x.max() - x.min()))
    data.loc[data["time"] <= 10., "Pdiff"] = np.nan
    data["Pmaxloc"] = data.apply(lambda row: data[(data["time_bin"] == row["time_bin"]) & (data["Pdiff"] <= Pthresh)]["Plith"].median() if pd.notna(row["Pdiff"]) and row["Pdiff"] <= Pthresh else np.nan, axis=1)

    # Step 1: Mark non-NaN Pmaxloc values
    data['Pmaxloc_not_nan'] = data['Pmaxloc'].notna()

    # Step 2: Create time_bin_consecutive
    data['time_bin_consecutive'] = np.where(
        data['Pmaxloc_not_nan'],
        (data['Pmaxloc_not_nan'] & (data['time_bin'] != data['time_bin'].shift())).cumsum(),
        np.nan
    )
    #if time_bin_consecutive is numeric and goes from value to value+1, use value instead of value+1
    for i in range(0, len(data)-1):
        if data["time_bin_consecutive"].iloc[i+1] > data["time_bin_consecutive"].iloc[i]:
            if abs(data["Pmaxloc"].iloc[i] - data["Pmaxloc"].iloc[i+1]) < Pthresh:
                data["time_bin_consecutive"].iloc[i+1] = data["time_bin_consecutive"].iloc[i]
    
    # see if the particle is rising, subducting or truly stagnant
    lastPdata = data["Plith"].iloc[-1] - data["Plith"].iloc[-10]
    if lastPdata <= - data["Plith"].max()/100.: 
        data["direction"] = "rising"
    elif lastPdata > data["Plith"].max()/100.:
        data["direction"] = "subducting"
    else:
        data["direction"] = "stagnant"
    
    # calculate difference in Plith between first value of one time_bin_consecutive_new and the last value of the next time_bin_consecutive_new
    data['Pdiff_new'] = data.groupby('time_bin_consecutive')['Plith'].transform(lambda x: x.iloc[-1] - x.iloc[0])
    data['Pmaxloc_new'] = data.groupby('time_bin_consecutive')['Pmaxloc'].transform("mean")
    return data




def main():
    parser = argparse.ArgumentParser(description= 'Script that gets some models and gives the kinematic indicators')
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


    for p in tqdm(range(0, len(part))):
        id = part["id"].iloc[p]
        data = load_data(id, txt_loc)

        Pthresh = data["Plith"].max()/20

        # Bin data by time intervals of 10 Myr
        bins = np.arange(0, 51, 5)
        data["time_bin"] = pd.cut(data["time"], bins)

        data = compute_stagnant_maxima(data, Pthresh)

        #save to particles dataframe
        particles["id"].iloc[p] = id
        particles["lithology"].iloc[p] = part["lithology"].iloc[p]
        particles["direction"].iloc[p] = data["direction"].iloc[0]
        particles["timePmax"].iloc[p] = data["time"][data["Plith"].idxmax()]
        particles["Pmax"].iloc[p] = data["Plith"].max()
        particles["Tmax"].iloc[p] = data["T"][data["Plith"].idxmax()]

        for k in range(len(data["Pmaxloc_new"].dropna().unique())): 
            particles[f"P{k+1}"].iloc[p] = data["Pmaxloc_new"].dropna().unique()[k]
            particles[f"timeP{k+1}"].iloc[p] = data[data["time_bin_consecutive"] == data["time_bin_consecutive"].dropna().unique()[k]]["time"].median()
            particles[f"T{k+1}"].iloc[p] = data[data["time_bin_consecutive"] == data["time_bin_consecutive"].dropna().unique()[k]]["T"].mean()

        if p % 10 == 0:
            plt.plot(data["time"], data["Plith"], color="grey")
            plt.plot(data["time"], data["Pmaxloc_new"], color="blue")
            plt.xlabel("Time (Myr)")
            plt.ylabel("Pressure (GPa)")
            plt.title(f"Particle {id} - Lithology: {part['lithology'].iloc[p]} - Direction: {data['direction'].iloc[0]}")
            plt.ylim(0, 2.0)
            plt.savefig(f"{sloc}/stagnant_maxima_{id}.png")
            plt.close()

    # particles = particles.astype({col: 'float' for col in particles.columns if col != 'lithology' or col != 'direction'})
    particles.to_csv(f"{txt_loc}/stagnant_maxima.txt", sep="\t", index=False, header=True, float_format='%.2f')




if __name__ == "__main__":
    main()