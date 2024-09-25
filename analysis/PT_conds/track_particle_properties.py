#!/usr/bin/env python3

# script that finds a particle through its index and collects the position and p-T path in time

import numpy as np
import os
import sys
import pandas as pd
import argparse
from pathlib import Path
path = str(Path(Path(__file__).parent.absolute()).parent.absolute())
sys.path.insert(0, path)
from libraries.functions import *
from libraries.particles import *
import matplotlib.pyplot as plt
import json
import seaborn as sns
from tqdm import tqdm
import warnings
warnings.filterwarnings("ignore") 


def main():

    parser = argparse.ArgumentParser(description= 'Script that gets n models and the time at the end of computation and gives the temperature and viscosity plots for all time steps')
    parser.add_argument('json_file', help='json file with models and end time to plot T and eta field')  
    args = parser.parse_args()

    csvs_loc =  '/home/vturino/PhD/projects/exhumation/gz_outputs/'
    models_loc =  '/home/vturino/PhD/projects/exhumation/raw_outputs/'
    json_loc = '/home/vturino/PhD/projects/exhumation/pyInput/'
 

    with open(f"{json_loc}{args.json_file}") as json_file:
            configs = json.load(json_file)

    m = configs['models'][0]
    plt_loc = f'/home/vturino/PhD/projects/exhumation/plots/single_models/{m}/'
    txt_loc = f'/home/vturino/PhD/projects/exhumation/plots/single_models/{m}/txt_files'
    pt_files = f'/home/vturino/PhD/projects/exhumation/plots/single_models/{m}/txt_files/tracked_particles'
    if not os.path.exists(pt_files):
        os.mkdir(pt_files)


    ap = pd.read_csv(f"{txt_loc}/all_particles_indexes.csv", sep='\t')
    time_array = np.zeros((len(os.listdir(f"{csvs_loc}{m}/fields")),2))   
    stat = pd.read_csv(f"{models_loc}{m}/statistics",skiprows=configs['head_lines'],sep='\s+',header=None)
    time_array = grab_dimTime_particles(f"{csvs_loc}{m}/fields", stat, time_array)
    ts = int(len(time_array)/2)

    # part0 =  pd.read_parquet(f"{csvs_loc}{m}/particles/full.0.gzip")
    # part0 = part0.reset_index()
    # merged = pd.merge(ap, part0, on=['initial position:0', 'initial position:1'], how='inner')
   

    # part1 =  pd.read_parquet(f"{csvs_loc}{m}/particles/full.1.gzip")
    # part1 = part1.reset_index()
    # merged1 = pd.merge(merged, part1, on=['index'], how='inner')
    # print(merged1)

    # missing_indexes = part0[~part0["index"].isin(ap['index'])]
    # if not missing_indexes.empty:
    #     print("Some indexes in part0 are not in ap:", missing_indexes)
    # else:
    #     print("All indexes in part0 are in ap.")


    # for t in range(2):
    #     particles = pd.read_parquet(f"{csvs_loc}{m}/particles/full.{t}.gzip")
    #     particles['terrain'] = np.where((particles['oc'] < particles['sed']), 'sed', np.where((particles['oc'] > particles['sed']), 'oc', 'none'))
    #     particles = particles.reset_index()
    #     tmp = particles[particles["index"].isin(part0['index'])].reset_index()
    #     print(tmp)
        # if tmp["index"].iloc[0] == part0["index"].iloc[0]:
        #     if tmp["initial position:0"].iloc[0] == part0["initial position:0"].iloc[0] and tmp["initial position:1"].iloc[0] == part0["initial position:1"].iloc[0]:
        #         print("For t =", t, "tmp['index'].iloc[0] has the same 'initial position:0' and 'initial position:1'")
        #     else:
        #         print("For t =", t, "tmp['index'].iloc[0] does not have the same 'initial position:0' and 'initial position:1'")
        # else:
        #     print("For t =", t, "tmp['index'].iloc[0] is not the same as part0['index'].iloc[0']")

    # for i in range(len(part0)):
    #         pt = open(f"{pt_files}/pt_part_{i}.txt", "w+")
    #         pt.write("ts x y xi yi P T terrain\n")
    # pt.close()

    # for t in tqdm(range(ts)):
    #     particles = pd.read_parquet(f"{csvs_loc}{m}/particles/full.{t}.gzip")
    #     particles = particles.reset_index()
    #     
    #     for i in range(len(tmp)):
    #         pt = open(f"{pt_files}/pt_part_{i}.txt", "a+")
    #         pt.write("%.0f %.3f %.3f %.3f %.3f %.3f %.3f %s \n" % (t, tmp["initial position:0"].iloc[i], tmp["initial position:1"].iloc[i], tmp["Points:0"].iloc[i], tmp["Points:1"].iloc[i], tmp["p"].iloc[i], tmp["T"].iloc[i], tmp["terrain"].iloc[i]))
    #     pt.close()
            



if __name__ == "__main__":
    main()