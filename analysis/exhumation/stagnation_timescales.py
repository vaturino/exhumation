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

    # csvs_loc =  '/home/vturino/PhD/projects/exhumation/gz_outputs/'
    # models_loc =  '/home/vturino/PhD/projects/exhumation/raw_outputs/'

    models_loc = '/home/vturino/Vale_nas/exhumation/raw_outputs/'
    csvs_loc = '/home/vturino/Vale_nas/exhumation/gz_outputs/'
    json_loc = '/home/vturino/PhD/projects/exhumation/pyInput/'
    rocks_loc = '/home/vturino/PhD/projects/exhumation/rock_record/'

    rocks = pd.read_excel(f"{rocks_loc}rocks_agard2018.xlsx")

    # Read the json file
    with open(f"{json_loc}{args.json_file}") as json_file:
            configs = json.load(json_file)
            setting = 'none'
    compositions = configs['compositions']
    m = configs["models"][0]

    # Create the folders to save the plots
    plot_loc = f"/home/vturino/PhD/projects/exhumation/plots/single_models/{configs['models'][0]}"
    txt_loc = f'{plot_loc}/txt_files'
    if not os.path.exists(txt_loc):
        os.mkdir(txt_loc)

    sloc = f"{plot_loc}/stagnant"
    if not os.path.exists(sloc):
        os.mkdir(sloc)

    stagloc = f"{txt_loc}/stagnant"
    if not os.path.exists(stagloc):
        os.mkdir(stagloc)

    pt_files = f'{txt_loc}/PT'
    npa = len(os.listdir(pt_files))

     #convergence rate
    cr = pd.read_csv(f"{txt_loc}/2D_v.txt", sep="\s+")
    cr["conv_rate"].iloc[0]= np.nan

    init = pd.read_csv(f"{txt_loc}/stagnant_particles.txt", sep="\s+")
    stag = pd.DataFrame(columns=["id", "maxP", "maxT", "lithology", "tin", "tmax", "tfin", "burial", "stag", "maxdepth", "stagdepth", "vbur", "vstag"], index=range(len(init)))

    print("Number of stagnant particles = ", len(init)) 
    
    for p in tqdm(range(len(init))):
    # for p in range(10):
        id = init["id"].iloc[p]
        part = pd.read_csv(f"{pt_files}/pt_part_{id}.txt", sep="\s+")
        part["Plith"] = (part["depth"].max()- part["depth"])*1e3*9.81*3300/1e9
        # print(part["Plith"])

        stag["id"].iloc[p] = id
        stag["maxP"].iloc[p] = init["maxPP"].iloc[p]
        stag["maxT"].iloc[p] = init["maxPT"].iloc[p]
        stag["lithology"].iloc[p] = init["lithology"].iloc[p]
        stag["tmax"].iloc[p] = init["tmax"].iloc[p]
        stag["maxdepth"] = part["depth"].max() - part["depth"].min()
        stag["stagdepth"] = 0.2*stag["maxdepth"]


        # stag["tin"].iloc[p] = part.loc[part["Plith"] >= 0.1, "time"].min()/2.
        stag["tfin"].iloc[p] = part["time"].iloc[-1]/2.
        for i in range(len(part)):
            if part.depth.iloc[0] - part.depth.iat[i] >= 2.:
                stag["tin"].iloc[p] = part["time"].iat[i]/2.
                break
    
    stag["burial"] = stag["tmax"] - stag["tin"]    
    stag["stag"] = stag["tfin"] - stag["tmax"]
    stag["vbur"] = stag["maxdepth"]*1.e5/(stag["burial"]*1.e6)
    stag["vstag"] = stag["stagdepth"]*1.e5/(stag["stag"]*1.e6)
            
    stag.to_csv(f"{txt_loc}/timing_stagnant_particles.txt", sep=" ", index=False)  

    f1, a1 = plt.subplots(2, 2, figsize=(10, 10))      
    sns.scatterplot(data=stag, x="maxT", y="maxP", hue="lithology", size = "vstag", ax=a1[0,0], zorder=10)
    a1[0,0].set_ylabel("Pressure (GPa)")
    a1[0,0].set_xlabel("T ($^\circ$C)")
    a1[0,0].set_title("Peak pressure")
    sns.kdeplot(data=rocks, x="T", y="P", ax=a1[0,0], fill = True, cbar = False, alpha = .7, zorder=1, color='grey')

    sns.histplot(x = "tmax", bins = 20, hue = "lithology", element="step", data=stag, ax=a1[0,1])
    a1[0,1].set_title("Time at peak pressure")
    a1[0,1].set_xlabel("Time (Ma)")
    a1[0,1].set_ylabel("Number of particles")
    ax1 = a1[0,1].twinx()
    ax1.plot(cr["time"]/1e6, cr["conv_rate"], color="grey", linewidth=2)
    ax1.set_ylabel("Convergence rate (cm/yr)", color="grey", fontweight='bold')
    a1[0,1].patch.set_visible(False) 
    a1[0,1].set_zorder(10)   
    

    sns.histplot(x = "vbur", bins = 20, hue = "lithology", element="step", data=stag, ax=a1[1,0])
    a1[1,0].set_title("Burial rate")
    a1[1,0].set_xlabel("Rate (cm/yr)")
    a1[1,0].set_ylabel("Number of particles")
    # a1[1,0].set_yscale('log')

    sns.histplot(x = "vstag", bins = 20, hue = "lithology", element="step", data=stag, ax=a1[1,1])
    a1[1,1].set_title("Stagnation rate")
    a1[1,1].set_xlabel("Rate (cm/yr)")
    a1[1,1].set_ylabel("Number of particles")
    # a1[1,1].set_yscale('log')

    f1.tight_layout()
    plt.savefig(f"{plot_loc}/timing_stagnant_particles.png", dpi = 1000)
    plt.close()
         



if __name__ == '__main__':
    main()