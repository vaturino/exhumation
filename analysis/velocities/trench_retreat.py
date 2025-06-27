#! /usr/bin/python3
import numpy as np
import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from numpy.core.fromnumeric import size
from scipy.interpolate import griddata
from matplotlib.gridspec import GridSpec
import sys, os, subprocess
import json as json
from tqdm import tqdm
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'
import argparse
from pathlib import Path
import seaborn as sns
import matplotlib.tri as tri
path = str(Path(Path(__file__).parent.absolute()).parent.absolute())
sys.path.insert(0, path)
from libraries.functions import *
from libraries.particles import *
import plotly.graph_objects as go
from itertools import takewhile


def get_trench_index(p, threshold = 0.3e7):
    if {"oc","sed"}.issubset(p.columns):
        sumcomp = p["oc"] + p["sed"]
        tr =  p.loc[(p['Points:0']> threshold) & (sumcomp > 0.3) & (p["Points:1"] >= p["Points:1"].max() - 5.e3),'Points:0'].idxmax()
    elif {"oc","serp"}.issubset(p.columns):
        sumcomp = p["oc"] + p["serp"]
        tr =  p.loc[(p['Points:0']> threshold) & (sumcomp > 0.3) & (p["Points:1"] >= p["Points:1"].max() - 5.e3),'Points:0'].idxmax()
    else:
        tr =  p.loc[(p['Points:0']> threshold) & (p['oc'] > 0.3) & (p["Points:1"] >= p["Points:1"].max() - 5.e3),'Points:0'].idxmax()
    return tr

def main():
    parser = argparse.ArgumentParser(description= 'Script that gets n models and the time at the end of computation and gives the temperature and viscosity plots for all time steps')
    parser.add_argument('json_file', help='json file with models and end time to plot T and eta field')  
    args = parser.parse_args()

    # csvs_loc =  '/home/vturino/PhD/projects/exhumation/gz_outputs/'
    # models_loc =  '/home/vturino/PhD/projects/exhumation/raw_outputs/'
    json_loc = '/home/vturino/PhD/projects/exhumation/pyInput/'

    models_loc = '/home/vturino/Vale_nas/exhumation/raw_outputs/'
    csvs_loc = '/home/vturino/Vale_nas/exhumation/gz_outputs/'

    with open(f"{json_loc}{args.json_file}") as json_file:
            configs = json.load(json_file)

    # compositions = configs['compositions']
    # cutoff = configs['cutoff']

    file_count = 0


    for ind_m, m in tqdm(enumerate(configs['models'])):    
        time_array = np.zeros((len(os.listdir(f"{csvs_loc}{m}/fields")),2)) 
        stat = pd.read_csv(f"{models_loc}{m}/statistics",skiprows=configs['head_lines'],sep='\s+',header=None)
        time_array = grab_dimTime_fields(f"{csvs_loc}{m}/fields", stat, time_array, configs['head_lines']-1)
        plot_loc_mod = f"/home/vturino/PhD/projects/exhumation/plots/single_models/{m}"
        if not os.path.exists(plot_loc_mod):
            os.mkdir(plot_loc_mod)
        plot_loc = f"{plot_loc_mod}/kinematics/"
        if not os.path.exists(plot_loc):
            os.mkdir(plot_loc)
        else:
            # Count the files in the fields_loc directory
            file_count = len(os.listdir(plot_loc))

        
        txt_loc = f"{plot_loc_mod}/txt_files/"
        if not os.path.exists(txt_loc):
            os.mkdir(txt_loc)
        else:
            # Count the files in the fields_loc directory
            file_count = len(os.listdir(txt_loc))

        
        
        ymax_plot = 902.e3

        
        plotname = "trench_motion.png"
        trench = pd.DataFrame(columns=["time", "trench_x", "trench_y"], index=range(len(time_array)))
        
        for t in tqdm(range(0, len(time_array))):
        # for t in tqdm(range(0,5)):
            data = pd.read_parquet(f"{csvs_loc}{m}/fields/full.{int(t)}.gzip") 

            pts = get_points_with_y_in(data, 15.e3, 2.e3, ymax=900.e3)
            trench_index = get_trench_index(pts, threshold=0.13e7)
            trench.loc[t, "time"] = time_array[t, 1]/1.e6
            trench.loc[t, "trench_x"] = data.loc[trench_index, "Points:0"]/1.e3
            trench.loc[t, "trench_y"] = (ymax_plot - data.loc[trench_index, "Points:1"])/1.e3
        
        trench.to_csv(f"{txt_loc}/trench_pos.csv", index=False)
             
        plt.plot(trench["time"], trench["trench_x"]-trench["trench_x"].iloc[0], label="Trench x position", color="blue", marker="o", markersize=2, linewidth=1)
        plt.axhline(0, color='black', linestyle='--', linewidth=0.5)
        plt.xlabel("Time (Ma)")
        plt.ylabel("Trench motion relative to initial (km)")
        plt.xlim(0, 50)
        # plt.ylim(0, 5400)
        plt.title(f"Trench position")
        plt.savefig(f"{plot_loc}{plotname}", dpi=500)
        plt.close("all")





if __name__ == "__main__":
    main()