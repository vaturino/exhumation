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

    eloc = f"{plot_loc}/exhumed"
    if not os.path.exists(eloc):
        os.mkdir(eloc)

    exhloc = f"{txt_loc}/exhumed"
    if not os.path.exists(exhloc):
        os.mkdir(exhloc)

    pt_files = f'{txt_loc}/PT'
    npa = len(os.listdir(pt_files))
    print("Total number of particles = ", npa)

    init = pd.read_csv(f"{txt_loc}/exhumed_particles.txt", sep="\s+")
    exh = pd.DataFrame(columns=["id", "maxP", "maxT", "lithology", "tin", "tmax", "tfin", "burial", "exhum", "maxdepth", "exdepth", "vbur", "vexh"], index=range(len(init)))
    
    for p in tqdm(range(len(init))):
    # for p in range(10):
        id = init["id"].iloc[p]
        part = pd.read_csv(f"{pt_files}/pt_part_{id}.txt", sep="\s+")
        part["Plith"] = (part["depth"].max()- part["depth"])*1e3*9.81*3300/1e9

        exh["id"].iloc[p] = id
        exh["maxP"].iloc[p] = init["maxPP"].iloc[p]
        exh["maxT"].iloc[p] = init["maxPT"].iloc[p]
        exh["lithology"].iloc[p] = init["lithology"].iloc[p]
        exh["tmax"].iloc[p] = init["tmax"].iloc[p]
        exh["maxdepth"].iloc[p] = (part["depth"].max() - part["depth"].min())
        exh["exdepth"].iloc[p] = 0.2*exh["maxdepth"].iloc[p]


        # exh["tin"].iloc[p] = part.loc[part["Plith"] >= 0.1, "time"].min()/2.
        

        for i in range(len(part)):
            if part.depth.iloc[0] - part.depth.iat[i] >= 2.:
                exh["tin"].iloc[p] = part["time"].iat[i]/2.
                break
        for i in range(len(part)): 
            idxmin = part.depth.idxmin()
            if i > idxmin:
                if part.depth.iat[i] - part.depth.min() >= exh["exdepth"].iloc[p]:
                    exh["tfin"].iloc[p] = part["time"].iat[i]/2.
                    break
    
    
    exh["burial"] = exh["tmax"] - exh["tin"]    
    exh["exhum"] = exh["tfin"] - exh["tmax"]
    exh["vbur"] = exh["maxdepth"]*1e5/(exh["burial"]*1e6)
    exh["vexh"] = exh["exdepth"]*1e5/(exh["exhum"]*1e6)
            
    exh.to_csv(f"{txt_loc}/timing_exhumed_particles.txt", sep=" ", index=False)  

    f1, a1 = plt.subplots(2, 2, figsize=(10, 10))      
    sns.scatterplot(data=exh, x="maxT", y="maxP", hue="lithology", size = "vexh",ax=a1[0,0])
    a1[0,0].set_ylabel("Pressure (GPa)")
    a1[0,0].set_xlabel("T ($^\circ$C)")
    a1[0,0].set_title("Peak pressure")

    sns.histplot(x = "tmax", bins = 20, hue = "lithology", element="step", data=exh, ax=a1[0,1])
    a1[0,1].set_title("Time at peak pressure")
    a1[0,1].set_xlabel("Time (Ma)")
    a1[0,1].set_ylabel("Number of particles")
    

    sns.histplot(x = "vbur", bins = 20, hue = "lithology", element="step", data=exh, ax=a1[1,0])
    a1[1,0].set_title("Burial rate")
    a1[1,0].set_xlabel("Rate (cm/yr)")
    a1[1,0].set_ylabel("Number of particles")

    sns.histplot(x = "vexh", bins = 20, hue = "lithology", element="step", data=exh, ax=a1[1,1])
    a1[1,1].set_title("Exhumation rate")
    a1[1,1].set_xlabel("Rate (cm/yr)")
    a1[1,1].set_ylabel("Number of particles")

    f1.tight_layout()
    plt.savefig(f"{plot_loc}/timing_exhumed_particles.png", dpi = 1000)
    plt.close()
         



if __name__ == '__main__':
    main()