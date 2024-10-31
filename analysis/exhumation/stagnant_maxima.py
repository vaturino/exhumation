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

    for p in tqdm(0, range(len(part)), 10):
        id = part[p]
        data = pd.read_csv(f"{txt_loc}/PT/pt_part_{id}.txt", sep="\s+")
        data["Plith"] = (901.e3 - data["y"])*3100*9.81/1e9
        data["time"] = data["time"]/2
        data["ts"] = data["time"]*2


        tin = 0
        for i in range(len(data)):
            if (data.y.iloc[0] - data.y.iloc[i]) > 2.:
                tin = data.time.iloc[i]
                break

        # Bin data by time intervals of 10 Myr
        bins = np.arange(0, 51, 10)
        data["time_bin"] = pd.cut(data["time"], bins)

        data["Pdiff"] = data.groupby("time_bin")["Plith"].transform(lambda x: abs(x.max() - x.min()))
        data.loc[data["time"] <= tin, "Pdiff"] = np.nan
        data[data["Pdiff"] <= data["Plith"].max()/10]
        data["Pmaxloc"] = data.apply(lambda row: data[(data["time_bin"] == row["time_bin"]) & (data["Pdiff"] <= 0.15)]["Plith"].median() if pd.notna(row["Pdiff"]) and row["Pdiff"] <= 0.15 else np.nan, axis=1)

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
                if abs(data["Pmaxloc"].iloc[i] - data["Pmaxloc"].iloc[i+1]) < data["Plith"].max()/10:
                    data["time_bin_consecutive"].iloc[i+1] = data["time_bin_consecutive"].iloc[i]
        
        # calculate difference in Plith between first value of one time_bin_consecutive_new and the last value of the next time_bin_consecutive_new
        data['Pdiff_new'] = data.groupby('time_bin_consecutive')['Plith'].transform(lambda x: x.iloc[-1] - x.iloc[0])
        data['Pmaxloc_new'] = data.groupby('time_bin_consecutive')['Pmaxloc'].transform("mean")

        #save to particles dataframe
        particles["id"].iloc[p] = id
        particles["timePmax"].iloc[p] = data["time"][data["Plith"].idxmax()]
        particles["Pmax"].iloc[p] = data["Plith"].max()

        for k in range(len(data["Pmaxloc_new"].dropna().unique())): 
            particles[f"P{k+1}"].iloc[p] = data["Pmaxloc_new"].dropna().unique()[k]
            particles[f"timeP{k+1}"].iloc[p] = data[data["time_bin_consecutive"] == data["time_bin_consecutive"].dropna().unique()[k]]["time"].median()
        
        plt.plot(data["time"], data["Plith"])
        plt.plot(data["time"], data["Pmaxloc_new"])
        plt.savefig(f"{sloc}/particle_{id}.png", dpi=300)
        plt.close()



    particles = particles.astype(float)
    particles["id"] = particles["id"].astype(int)
    particles.to_csv(f"{txt_loc}/stagnant_maxima.txt", sep="\t", index=False, header=True, float_format='%.2f')




if __name__ == "__main__":
    main()