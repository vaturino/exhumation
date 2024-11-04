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

def find_first_stable_gradient_threshold(depth, tmax, threshold=0.1):
    # Calculate the gradient of the depth array
    gradient = np.gradient(depth)

    # Find the first index where the gradient stabilizes (>= -threshold) and remains stable until the end
    for i in range(len(gradient)):
        if np.all(gradient[i:] >= -threshold) and np.all(gradient[i:] <= threshold) and (i/2.) > tmax:
            return i
    
    return None  # Return None if no such index is found

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
    rocks = rocks[rocks["AREA"] != "Metamorphic Soles"]

    colors_tin = {
        "sed": "midnightblue",
        "oc": "#733A11",
        "ecl": "#003300",
        "serp": "#3b0000"
    }

    colors_tmax = {
        "sed": "mediumblue",
        "oc": "#B06D1A",
        "ecl": "#45701C",
        "serp": "brown"
    }

    colors_tfin = {
        "sed": "cornflowerblue",
        "oc": "#E3B64F",
        "ecl": "#A0C93D",
        "serp": "lightsalmon"
    }

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
    subd = pd.read_csv(f"{txt_loc}/subducted_particles.txt", sep="\s+")
    subd = subd.dropna()

    stag = pd.DataFrame(columns=["id", "maxP", "maxT", "lithology", "tin", "tmax", "tfin", "burial", "stag", "maxdepth", "stagdepth", "vbur", "vstag", "Pin", "Pstag"], index=range(len(init)))

    print("Number of stagnant particles = ", len(init)) 

    ymax = 900.
    
    for p in tqdm(range(len(init))):
    # for p in range(10):
        id = init["id"].iloc[p]
        part = pd.read_csv(f"{pt_files}/pt_part_{id}.txt", sep="\s+")

        stag["id"].iloc[p] = id
        stag["maxP"].iloc[p] = init["maxPP"].iloc[p]
        stag["maxT"].iloc[p] = init["maxPT"].iloc[p]
        stag["lithology"].iloc[p] = init["lithology"].iloc[p]
        stag["tmax"].iloc[p] = init["tmax"].iloc[p]
        stag["maxdepth"].iloc[p] = (ymax - part["depth"].min())
        stag["stagdepth"] = 0.2*stag["maxdepth"]

        idx = 0
        ide = 0

        finidx = find_first_stable_gradient_threshold(part["depth"], stag["tmax"].iloc[p], threshold=0.1)
        if finidx is not None:
            stag["tfin"].iloc[p] = part["time"].iat[finidx]/2.
            stag["Pstag"].iloc[p] = part["Plith"].iat[finidx]
            ide = finidx
        else:
            stag["tfin"].iloc[p] = part["time"].iat[-2]/2.
            stag["Pstag"].iloc[p] = part["Plith"].iat[-2]
            ide = len(part)

        for i in range(len(part)):
            if (part.depth.iloc[0] - part.depth.iloc[i]) >= 2.:
                stag["tin"].iloc[p] = part["time"].iat[i]/2.
                stag["Pin"].iloc[p] = part["Plith"].iat[i]
                idx = i
                break


        # vel_bur = part.loc[idx:part.depth.idxmin(), ["vx", "vy"]].apply(lambda row: math.sqrt(row["vx"]**2 + row["vy"]**2)*np.sin(np.deg2rad(45)), axis=1).values
        # vel_stag = part.loc[part.depth.idxmin():, ["vx", "vy"]].apply(lambda row: math.sqrt(row["vx"]**2 + row["vy"]**2)*np.sin(np.deg2rad(45)), axis=1).values
        # vel_bur = part.loc[idx:part.depth.idxmin(), ["vy"]].apply(lambda row: math.sqrt(row["vy"]**2), axis=1).values
        # vel_stag = part.loc[ide:, ["vy"]].apply(lambda row: math.sqrt(row["vy"]**2), axis=1).values

        # stag["vbur"].iloc[p] = vel_bur.mean()
        # stag["vstag"].iloc[p] = vel_stag.mean()
   
    
    stag["burial"] = stag["tmax"] - stag["tin"]   
    stag["stag"] = stag["tfin"] - stag["tmax"]
    stag["vbur"] = stag["maxdepth"]*1.e5/(stag["burial"]*1.e6)
    stag["vstag"] = stag["stagdepth"]*1.e5/(stag["stag"]*1.e6)
            
    stag.to_csv(f"{txt_loc}/timing_stagnant_particles.txt", sep=" ", index=False)  

    sns.set_palette("colorblind")

    f1, a1 = plt.subplots(2, 2, figsize=(10, 10))     

    for s in range(0, len(subd), 50):
        id = subd["id"].iloc[s].astype(int)
        spart = pd.read_csv(f"{pt_files}/pt_part_{id}.txt", sep="\s+")
        a1[0,1].plot((spart["T"]), spart["Plith"], color = 'darkslateblue', alpha = 0.5, linewidth = 0.5, zorder = 20)
            # a1[0,1].scatter(spart["T"].iat[-1]-273.15, spart["Plith"].iat[-1], color = 'darkslateblue', s = 10, zorder = 20)
    a1[0,1].set_ylabel("Pressure (GPa)")
    a1[0,1].set_xlabel("Temperature (C)")
    a1[0,1].set_title("Subducted particles")
    sns.kdeplot(data=rocks, x="T", y="P", ax=a1[0,1], fill = True, cbar = False, alpha = .7, zorder=1, color='grey')
    a1[0,1].set_ylim(0, 3.5)
    a1[0,1].set_xlim(0, 900)
    # # Bring the line to the front
    # for artist in a1[0,1].get_children():
    #     if isinstance(artist, plt.Line2D):  # Ensure only lines are affected
    #         artist.set_zorder(0)  # Move histogram lines behind


    sns.scatterplot(data=stag, 
                    x="maxT", 
                    y="maxP", 
                    hue="lithology", 
                    hue_order=stag["lithology"].value_counts(ascending=True).index,  
                    ax=a1[0,0], 
                    zorder = 10,
                    alpha=1)
    a1[0,0].set_ylabel("Pressure (GPa)")
    a1[0,0].set_xlabel("T ($^\circ$C)")
    a1[0,0].set_title("Peak pressure")
    a1[0,0].set_xlim(0, 900)
    a1[0,0].set_ylim(0, 3.5)
    sns.kdeplot(data=rocks, x="T", y="P", ax=a1[0,0], fill = True, cbar = False, alpha = .7, zorder=2, color='grey')
    


    sns.histplot(x = "tmax", 
                 bins = 20, 
                 hue = "lithology", 
                 hue_order=stag["lithology"].value_counts(ascending=True).index, 
                 element="step", 
                 data=stag, 
                 ax=a1[1,0], 
                 alpha=1,
                 edgecolor='black',
                 linewidth=1,
                 zorder=1)
    a1[1,0].set_title("Time at peak pressure")
    a1[1,0].set_xlabel("Time (Ma)")
    a1[1,0].set_ylabel("Number of particles")
    ax1 = a1[1,0].twinx()
    ax1.plot(cr["time"]/1e6, cr["conv_rate"], color="grey", linewidth=2, zorder =10)
    ax1.set_ylabel("Convergence rate (cm/yr)", color="grey", fontweight='bold')
    a1[1,0].patch.set_visible(False) 
    a1[1,0].set_zorder(1) 
    ax1.set_zorder(10)

    # Bring the line to the front
    for artist in a1[1,0].get_children():
        if isinstance(artist, plt.Line2D):  # Ensure only lines are affected
            artist.set_zorder(0)  # Move histogram lines behind  
    

    sns.histplot(x="vbur", 
                 bins=20, 
                 hue="lithology", 
                 hue_order=stag["lithology"].value_counts(ascending=True).index, 
                 element="step", 
                 data=stag, 
                 ax=a1[1,1], 
                 alpha=1, 
                 edgecolor='black', 
                 linewidth=1)
    a1[1,1].set_title("Burial rate")
    a1[1,1].set_xlabel("Rate (cm/yr)")
    a1[1,1].set_ylabel("Number of particles")
    # a1[1,0].set_yscale('log')

    # sns.histplot(x = "vstag", 
    #              bins = 20, 
    #              hue = "lithology", 
    #              hue_order=stag["lithology"].value_counts(ascending=True).index, 
    #              element="step", 
    #              data=stag, 
    #              ax=a1[1,1], 
    #              alpha=1,
    #              edgecolor='black',
    #              linewidth=1)
    # a1[1,1].set_title("Stagnation rate")
    # a1[1,1].set_xlabel("Rate (cm/yr)")
    # a1[1,1].set_ylabel("Number of particles")
    # # a1[1,1].set_yscale('log')

    f1.tight_layout()
    plt.savefig(f"{plot_loc}/timing_stagnant_particles.png", dpi = 1000)
    plt.close()
         



if __name__ == '__main__':
    main()