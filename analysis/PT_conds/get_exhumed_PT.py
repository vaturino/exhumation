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

    csvs_loc =  '/home/vturino/PhD/projects/exhumation/gz_outputs/'
    models_loc =  '/home/vturino/PhD/projects/exhumation/raw_outputs/'
    json_loc = '/home/vturino/PhD/projects/exhumation/pyInput/'
    cmyr = 1e2*(60*60*24*365)


    with open(f"{json_loc}{args.json_file}") as json_file:
            configs = json.load(json_file)
    
    for ind_m, m in tqdm(enumerate(configs['models'])): 
        time_array = np.zeros((len(os.listdir(f"{csvs_loc}{m}/fields")),2))   
        stat = pd.read_csv(f"{models_loc}{m}/statistics",skiprows=configs['head_lines'],sep='\s+',header=None)
        time_array = grab_dimTime_particles(f"{csvs_loc}{m}/fields", stat, time_array, configs['head_lines']-1)

        plot_loc = f"/home/vturino/PhD/projects/exhumation/plots/single_models/{configs['models'][ind_m]}"
        txt_loc = f'{plot_loc}/txt_files'
        if not os.path.exists(txt_loc):
            os.mkdir(txt_loc)
        pt_loc = f'{txt_loc}/part_indexes'
        if not os.path.exists(pt_loc):
            os.mkdir(pt_loc)

        parts = f"{pt_loc}/pt_part_0_Myr_merged.txt"
        ts = int(len(os.listdir(f"{pt_loc}")))
        with open(parts,"r") as f:
            npart = len(f.readlines()) -1
    

        P = np.zeros((npart, ts))
        T = np.zeros((npart, ts))
        P[:]=np.nan
        T[:]=np.nan
        pt_files = f'{txt_loc}/PT'
        if not os.path.exists(pt_files):
            os.mkdir(pt_files)


        for i in range(npart):
            pt = open(f"{pt_files}/pt_part_{i}.txt", "w+")
            pt.write("ts x y P T depth vy ocean sediments time_at_trench\n")
        pt.close()

        conds = pd.DataFrame(columns=["position:0", "position:1", "p", "T", "index", "velocity:1","oc", "sed"])
        for t in tqdm(range(1,ts)):
        # for t in tqdm(range(12,13)):
            idx = pd.read_csv(f"{pt_loc}/pt_part_{t}_Myr_merged.txt", sep="\s+", usecols= ["index", "time_at_trench"])
            df = pd.read_parquet(f"{csvs_loc}{m}/particles/full.{t}.gzip", columns= ["Points:0", "Points:1", "p", "T", "Points:1", "velocity:1", "oc", "sed"])
            df=df.reset_index()
            conds = pd.merge(idx, df, on=["index"], how="inner")
            # print(conds)
            P[:, t] = conds["p"]/1.e9
            T[:, t] = conds["T"]
            for i in range(npart):
                pt = open(f"{pt_files}/pt_part_{i}.txt", "a+")
                pt.write("%.0f %.3f %.3f  %.3f %.3f %.3f %.3f %.3f %.3f %.3f\n" % (t, conds["Points:0"].iloc[i], conds["Points:1"].iloc[i], P[i, t], T[i, t], conds["Points:1"].iloc[i]/1.e3, conds["velocity:1"].iloc[i]*cmyr, conds["oc"].iloc[i], conds["sed"].iloc[i], conds["time_at_trench"].iloc[i]))
            pt.close()
     


if __name__ == "__main__":
    main()


