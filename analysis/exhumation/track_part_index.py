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
import argparse
from pathlib import Path
import seaborn as sns
import matplotlib.tri as tri
path = str(Path(Path(__file__).parent.absolute()).parent.absolute())
sys.path.insert(0, path)
from libraries.functions import *
from libraries.particles import *
import plotly.graph_objects as go

def main():

    parser = argparse.ArgumentParser(description= 'Script that gets n models and the time at the end of computation and gives the temperature and viscosity plots for all time steps')
    parser.add_argument('json_file', help='json file with models and end time to plot T and eta field')  
    args = parser.parse_args()

    csvs_loc =  '/home/vturino/PhD/projects/exhumation/gz_outputs/'
    models_loc =  '/home/vturino/PhD/projects/exhumation/raw_outputs/'
    json_loc = '/home/vturino/PhD/projects/exhumation/pyInput/'

    with open(f"{json_loc}{args.json_file}") as json_file:
            configs = json.load(json_file)


    for ind_m, m in tqdm(enumerate(configs['models'])):    
        time_array = np.zeros((len(os.listdir(f"{csvs_loc}{m}/fields")),2)) 
        stat = pd.read_csv(f"{models_loc}{m}/statistics",skiprows=configs['head_lines'],sep='\s+',header=None)
        time_array = grab_dimTime_fields(f"{csvs_loc}{m}/fields", stat, time_array, configs['head_lines']-1)
        plot_loc_mod = f"/home/vturino/PhD/projects/exhumation/plots/single_models/{m}"
        if not os.path.exists(plot_loc_mod):
            os.mkdir(plot_loc_mod)
        plot_loc = f"{plot_loc_mod}/Density/"
        if not os.path.exists(plot_loc):
            os.mkdir(plot_loc)

        for t in tqdm(range(0, len(time_array), 2)):
        # for t in tqdm(range(2,3)):
            part = pd.read_parquet(f"{csvs_loc}{m}/particles/full.{int(t/2)}.gzip")
            fil = part["id"] == 205678
            # print("initial position = ", part[fil]["initial position:0"].to_numpy()/1e3, part[fil]["initial position:1"].to_numpy()/1e3)
            plt.plot(part[fil]["Points:0"].to_numpy()/1e3, part[fil]["Points:1"].to_numpy()/1e3, 'o', color='red')
        plt.savefig(f"{plot_loc_mod}/track_part.png")
        plt.close()

if __name__ == "__main__":
    main()