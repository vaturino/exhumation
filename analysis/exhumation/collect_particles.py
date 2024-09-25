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

    xmax = 5400.e3
    samples = 10
    zmax = 900.e3
    dlim =30.e3
    tr = 1.e-2


    with open(f"{json_loc}{args.json_file}") as json_file:
            configs = json.load(json_file)
    
    for ind_m, m in tqdm(enumerate(configs['models'])): 
        time_array = np.zeros((len(os.listdir(f"{csvs_loc}{m}/fields")),2))   
        stat = pd.read_csv(f"{models_loc}{m}/statistics",skiprows=configs['head_lines'],sep='\s+',header=None)
        time_array = grab_dimTime_particles(f"{csvs_loc}{m}/fields", stat, time_array, configs['head_lines']-1)

        plot_loc = f"../plots/single_models/{configs['models'][ind_m]}"
        txt_loc = f'{plot_loc}/txt_files'
        if not os.path.exists(txt_loc):
            os.mkdir(txt_loc)
        

        init = pd.read_parquet(f"{csvs_loc}{m}/particles/full.0.gzip")
        fildata = (init["oc"]>tr) | (init["sed"]>tr) | (init["opc"]>tr) | (init["ecl"]>tr)
        init = init[fildata]
        p = pd.DataFrame(columns=init.columns)
        p["time_at_trench"] = 0

        ts = int(len(time_array))

        # loop that, for every timestep, collects a number of particles "samples" that are at a distance of 50 km from the trench in the first "dlim" km of the domain
        for t in tqdm(range(0, ts-5)):
            full = load_dataset(loc_data=f"{csvs_loc}{m}/particles/full.{t}.gzip", data=init)
            data = pd.read_parquet(f"{csvs_loc}{m}/fields/full.{int(t)}.gzip")
            pts = get_points_with_y_in(data, 10.e3, 2.e3, ymax = 900.e3)
            trench= get_trench_position(pts,threshold = 0.13e7)
            incoming = trench -50.e3
            tmp = get_incoming_particles(full, 1., incoming, samples)
            tmp = tmp.assign(time_at_trench = t)
            p = pd.concat([p,tmp])
            # print(tmp)
        print("initial number of particles = ", len(p))

        # check that every particle is unique and identify the particles that are already in the first timestep
        p = p.drop_duplicates(subset = ["id"], keep='last')
        pi = pd.read_parquet(f"{csvs_loc}{m}/particles/full.{0}.gzip", columns=["id", "initial position:0", "initial position:1", "oc", "sed"]).reset_index()
        # pi = pi[pi['id'].isin(p['id'])]
        tp = pd.merge(p, pi, on="id", how="inner")
        npart = len(tp)
        print("total number of tracked particles = ", npart)
        # print(tp)
        plt.scatter(tp['Points:0'], tp['Points:1'], c=tp['oc'], cmap='viridis')
        plt.scatter(trench, 900.e3, c='red')
        plt.savefig(plot_loc + '/incoming_particles.png', dpi = 1000)
        plt.close()

        
        
        filename = f"{txt_loc}/particles_indexes.txt"
        pt = open(filename,"w+")
        pt.write("particle id \n")
        for i in range(0,npart):
            pt.write("%.0f %.3f\n" % (i, tp["id"].iloc[i]))
    
        pt.close()

       

        # # write the particles to a file
        # for t in tqdm(range(0,ts)):
        #     allp = pd.read_parquet(f"{csvs_loc}{m}/particles/full.{int(t)}.gzip", columns=["id", "initial position:0", "initial position:1", "oc", "sed"]).reset_index()
        #     merged = allp[allp['id'].isin(pi['id'])].sort_values(by=['id'])

        #     filename = f"{pt_loc}/pt_part_{t}_Myr_merged.txt"
        #     pt = open(filename,"w+")
        #     pt.write("particle id x y\n")
        #     for i in range(0,npart):
        #         pt.write("%.0f %.3f %.3f %.3f\n" % (i, merged["id"].iloc[i], merged["initial position:0"].iloc[i], merged["initial position:1"].iloc[i]))
            
        #     pt.close()


if __name__ == "__main__":
    main()


