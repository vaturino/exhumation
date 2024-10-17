#! /usr/bin/python3

import pandas as pd
import numpy as np
import json
from tqdm import tqdm
import matplotlib.pyplot as plt
import argparse
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import seaborn as sns
from matplotlib.cm import get_cmap

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
pd.options.mode.chained_assignment = None  # default='warn'

path = str(Path(Path(__file__).parent.absolute()).parent.absolute())
sys.path.insert(0, path)
from libraries.functions import *
from libraries.particles import *
from libraries.exhumation import *  



def get_trench_position_from_sed(p, threshold = 1.5e7):
    if {"sed"}.issubset(p.columns):
        tr =  p.loc[(p["sed"] > 0.3) & (p["Points:1"] >= p["Points:1"].max() - 10.e3),['Points:0', 'Points:1']].max()
    else:
        tr =  p.loc[(p['Points:0']> threshold) & (p['oc'] > 0.3),['Points:0', 'Points:1']].max()
    return tr


def main():
    parser = argparse.ArgumentParser(description='Script that gets some models and gives the kinematic indicators')
    parser.add_argument('json_file', help='json file with model name, time at the end of computation, folder where to save the plots, legend')
    args = parser.parse_args()

    models_loc = '/home/vturino/Vale_nas/exhumation/raw_outputs/'
    csvs_loc = '/home/vturino/Vale_nas/exhumation/gz_outputs/'
    json_loc = '/home/vturino/PhD/projects/exhumation/pyInput/'

    with open(f"{json_loc}{args.json_file}") as json_file:
        configs = json.load(json_file)
    compositions = configs['compositions']
    m = configs['models'][0]

    plot_loc = f"/home/vturino/PhD/projects/exhumation/plots/single_models/{m}"
    txt_loc = f'{plot_loc}/txt_files'
    os.makedirs(txt_loc, exist_ok=True)

    time_array = np.zeros((len(os.listdir(f"{csvs_loc}{m}/fields")),2)) 
    stat = pd.read_csv(f"{models_loc}{m}/statistics",skiprows=configs['head_lines'],sep='\s+',header=None)
    time_array = grab_dimTime_fields(f"{csvs_loc}{m}/fields", stat, time_array, configs['head_lines']-1)
    trench = pd.DataFrame(np.zeros((len(time_array),3)), columns=['time', 'x', 'y'])
    for t in tqdm(range(0, len(time_array))):
        p = pd.read_parquet(f"{csvs_loc}{m}/fields/full.{int(t)}.gzip") 
        trench.loc[t, 'time'] = time_array[t, 1]/1.e6
        trench.loc[t, 'x'], trench.loc[t, 'y'] = get_trench_position_from_sed(p)
        print(trench.loc[t, 'time'], trench.loc[t, 'x'], trench.loc[t, 'y'])
        
    trench.to_csv(f"{txt_loc}/trench_pos.txt", sep=' ', index=False)

    plt.scatter(trench['x']/1.e3, trench['y']/1.e3, c=trench['time']/1.e6, cmap='viridis')
    plt.colorbar()
    plt.xlabel('x [m]')
    plt.ylabel('y [m]')
    plt.savefig(f"{plot_loc}/trench_pos.png", dpi = 500)







if __name__ == "__main__":
    main()