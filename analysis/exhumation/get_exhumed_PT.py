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
path = str(Path(Path(__file__).parent.absolute()).parent.absolute())
sys.path.insert(0, path)
from libraries.functions import *
from libraries.particles import *


def main():

    parser = argparse.ArgumentParser(description= 'Script that gets some models and gives the kinematic indicators')
    parser.add_argument('json_file', help='json file with model name, time at the end of computation, folder where to save the plots, legend')
    args = parser.parse_args()

    # csvs_loc =  '/home/vturino/PhD/projects/exhumation/gz_outputs/'
    # models_loc =  '/home/vturino/PhD/projects/exhumation/raw_outputs/'

    models_loc = '/home/vturino/Vale_nas/exhumation/raw_outputs/'
    csvs_loc = '/home/vturino/Vale_nas/exhumation/gz_outputs/'

    json_loc = '/home/vturino/PhD/projects/exhumation/pyInput/'
    cmyr = 1e2*(60*60*24*365)
    tr = 1.e-2

    with open(f"{json_loc}{args.json_file}") as json_file:
            configs = json.load(json_file)

    compositions = configs['compositions']

    # Define the dataframe columns
    fixed_columns = ["position:0", "position:1", "p", "T", "id", "velocity:1"]
    columns = fixed_columns + compositions 

    
    for ind_m, m in tqdm(enumerate(configs['models'])): 
        time_array = np.zeros((len(os.listdir(f"{csvs_loc}{m}/particles")),2))   
        stat = pd.read_csv(f"{models_loc}{m}/statistics",skiprows=configs['head_lines'],sep='\s+',header=None)
        time_array = grab_dimTime_particles(f"{csvs_loc}{m}/particles", stat, time_array, configs['head_lines']-1)


        plot_loc = f"/home/vturino/PhD/projects/exhumation/plots/single_models/{configs['models'][ind_m]}"
        txt_loc = f'{plot_loc}/txt_files'
        if not os.path.exists(txt_loc):
            os.mkdir(txt_loc)

        # get number of tracked particles
        parts = f"{txt_loc}/particles_indexes.txt"
        ts = len(time_array)
        with open(parts,"r") as f:
            npart = len(f.readlines()) -1
    
        pt_files = f'{txt_loc}/PT'
        if not os.path.exists(pt_files):
            os.mkdir(pt_files)

        conds = pd.DataFrame(columns=columns)
        idx = pd.read_csv(f"{txt_loc}/particles_indexes.txt", sep="\s+")

        for i in range(npart):
            pt = open(f"{pt_files}/pt_part_{i}.txt", "w+")
            pt.write("id time x y P T depth vy " + " ".join(compositions) + "\n")
        pt.close()

        for t in tqdm(range(1,ts)):
            fcol = ["id", "position:0", "position:1", "p", "T", "position:1", "velocity:1"]
            df = pd.read_parquet(f"{csvs_loc}{m}/particles/full.{t}.gzip", columns= fcol + compositions)
            fil = (df[compositions]> tr).any(axis=1)
            df = df[fil]
            conds = df[df['id'].isin(idx['id'])].sort_values(by='id') 
      
            
        # Write filtered and sorted data to particle files
            for i in range(npart):
                pt = open(f"{pt_files}/pt_part_{i}.txt", "a+")
                pt.write("%.0f %.0f %.3f %.3f  %.3f %.3f %.3f %.3f " % (
                    conds["id"].iloc[i], t, conds["position:0"].iloc[i], conds["position:1"].iloc[i],
                    conds["p"].iloc[i] / 1e9, conds["T"].iloc[i], conds["position:1"].iloc[i] / 1.e3,
                    conds["velocity:1"].iloc[i] * cmyr
                ))
                
                # Dynamically write composition values
                for col in compositions:
                    if col == compositions[-1]:
                        pt.write("%.3f\n" % conds[col].iloc[i])
                    else:
                        pt.write("%.3f " % conds[col].iloc[i])
                
            pt.write("\n")
            pt.close()             
     


if __name__ == "__main__":
    main()


