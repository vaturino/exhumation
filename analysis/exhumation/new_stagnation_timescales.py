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
    parser = argparse.ArgumentParser(description='Script that gets some models and gives the kinematic indicators')
    parser.add_argument('json_file', help='json file with model name, time at the end of computation, folder where to save the plots, legend')
    parser.add_argument('--models_dir', help='Directory for raw models output', default='/home/vturino/Vale_nas/exhumation/raw_outputs/')
    parser.add_argument('--csvs_dir', help='Directory for CSV files', default='/home/vturino/Vale_nas/exhumation/gz_outputs/')
    parser.add_argument('--json_dir', help='Directory for JSON input', default='/home/vturino/PhD/projects/exhumation/pyInput/')
    parser.add_argument('--rocks_dir', help='Directory for rocks data', default='/home/vturino/PhD/projects/exhumation/rock_record/')
    args = parser.parse_args()

    # Use directories from arguments
    models_loc = args.models_dir
    csvs_loc = args.csvs_dir
    json_loc = args.json_dir
    rocks_loc = args.rocks_dir

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

    # Convergence rate
    cr = pd.read_csv(f"{txt_loc}/2D_v.txt", sep="\s+")
    cr["conv_rate"].iloc[0] = np.nan

    init = pd.read_csv(f"{txt_loc}/stagnant_particles.txt", sep="\s+")
    subd = pd.read_csv(f"{txt_loc}/subducted_particles.txt", sep="\s+")
    subd = subd.dropna()

    stag = pd.DataFrame(columns=["id", "Pm_kin", "Pm_dyn", "Pm_trans", "Tm_kin", "Tm_dyn", "Tm_trans", "lithology", "tin", "tm_kin", "tm_dyn", "tm_trans", "burial", "stag", "maxdepth", "stagdepth", "vbur", "vstag", "Pin", "Pstag"], index=range(len(init)))

    print("Number of stagnant particles = ", len(init)) 

    ymax = 900.
    
    for p in tqdm(range(len(init))):
        id = init["id"].iloc[p]
        part = pd.read_csv(f"{pt_files}/pt_part_{id}.txt", sep="\s+")

        if init["Pm_kin"].notna().iloc[p]:
            stag["Pm_kin"].iloc[p] = init["Pm_kin"].iloc[p]
            stag["Tm_kin"].iloc[p] = init["Tm_kin"].iloc[p]
            stag["tm_kin"].iloc[p] = init["tm_kin"].iloc[p]
        if init["Pm_dyn"].notna().iloc[p]:
            stag["Pm_dyn"].iloc[p] = init["Pm_dyn"].iloc[p]
            stag["Tm_dyn"].iloc[p] = init["Tm_dyn"].iloc[p]
            stag["tm_dyn"].iloc[p] = init["tm_dyn"].iloc[p]
        if init["Pm_trans"].notna().iloc[p]:
            stag["Pm_trans"].iloc[p] = init["Pm_trans"].iloc[p]
            stag["Tm_trans"].iloc[p] = init["Tm_trans"].iloc[p]
            stag["tm_trans"].iloc[p] = init["tm_trans"].iloc[p]

        stag["id"].iloc[p] = id
        stag["lithology"].iloc[p] = init["lithology"].iloc[p]
        stag["maxdepth"].iloc[p] = (ymax - part["depth"].min())
        stag["stagdepth"] = 0.25 * stag["maxdepth"]

        idx = 0
        ide = 0


        for i in range(len(part)):
            if (part.depth.iloc[0] - part.depth.iloc[i]) >= 2.:
                stag["tin"].iloc[p] = part["time"].iat[i]/2.
                stag["Pin"].iloc[p] = part["Plith"].iat[i]
                idx = i
                break

        # # Calculate rates and times
        # stag["burial"] = stag["tmax"] - stag["tin"]   
        # stag["stag"] = stag["tfin"] - stag["tmax"]
        # stag["vbur"] = stag["maxdepth"] * 1.e5 / (stag["burial"] * 1.e6)
        # stag["vstag"] = stag["stagdepth"] * 1.e5 / (stag["stag"] * 1.e6)
    
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

    print("plotted subducted particles")
        
    sns.kdeplot(data=rocks, x="T", y="P", ax=a1[0,0], fill = True, cbar = False, alpha = .7, zorder=1, color='grey')
    # for i in range(0, len(stag)):
    sns.scatterplot(data=stag, x="Tm_kin", y="Pm_kin", hue="lithology", ax=a1[0,0], palette=colors_tmax, legend=True)
    sns.scatterplot(data=stag, x="Tm_dyn", y="Pm_dyn", hue="lithology", ax=a1[0,0], palette=colors_tmax, legend=False)
    sns.scatterplot(data=stag, x="Tm_trans", y = "Pm_trans", hue="lithology", ax=a1[0,0], palette=colors_tmax, legend=False)
    a1[0,0].set_xlabel('Time (My)')
    a1[0,0].set_ylabel('Pressure (km)')
    a1[0,0].set_ylim(0, 3.5)
    a1[0,0].set_title('Stagnant Particle Trajectories')
    f1.tight_layout()

    # Figure 3: histogram of stagnation times tm_i
    sns.histplot(x = "tm_dyn", 
                 bins = 20, 
                 hue = "lithology", 
                 hue_order=stag["lithology"].value_counts(ascending=True).index, 
                 palette=colors_tmax,
                 element="step", 
                 data=stag, 
                 ax=a1[1,0], 
                 alpha=1,
                 edgecolor='black',
                 linewidth=1,
                 zorder=1)
    sns.histplot(x = "tm_kin",
                 bins = 20,
                 hue = "lithology",
                 hue_order=stag["lithology"].value_counts(ascending=True).index,
                 palette=colors_tmax,
                 element="step",
                 data=stag,
                 ax=a1[1,0],
                 alpha=1,
                 edgecolor='black',
                 linewidth=1,
                 zorder = 1)
    sns.histplot(x = "tm_trans",
                 bins = 20,
                 hue = "lithology",
                 hue_order=stag["lithology"].value_counts(ascending=True).index,
                 palette=colors_tmax,
                 element="step",
                 data=stag,
                 ax=a1[1,0],
                 alpha=1,
                 edgecolor='black',
                 linewidth=1,
                 zorder = 1)
    

    a1[1,0].set_title("Time at peak pressure")
    a1[1,0].set_xlabel("Time (Ma)")
    a1[1,0].set_ylabel("Number of particles")
    ax1 = a1[1,0].twinx()
    ax1.plot(cr["time"]/1e6, cr["conv_rate"], color="grey", linewidth=2, zorder =10)
    ax1.set_ylabel("Convergence rate (cm/yr)", color="grey", fontweight='bold')
    a1[1,0].patch.set_visible(False) 
    a1[1,0].set_zorder(1) 
    ax1.set_zorder(10)


    # Plot 4: pie chart of how many particles stagnate at dyn, trans, kin
    dyn = len(stag[stag["tm_dyn"].notna()])
    kin = len(stag[stag["tm_kin"].notna()])
    trans = len(stag[stag["tm_trans"].notna()])

    #pie chart
    a1[1,1].pie([dyn,kin, trans], labels = ['Dynamic\nslowdown', 'Kinematic\nslowdown', 'Transition\npoint'], colors = ['lightcoral', 'lightblue', 'blueviolet'], autopct='%1.1f%%', startangle=90, radius = 0.9)
    a1[1,1].text(-1.5, -1.1, "Number of particles stagnating during dynamic slowdown: " + str(dyn), fontsize=10)
    a1[1,1].text(-1.5, -1.2, "Number of particles stagnating during kinematic slowdown: " + str(kin), fontsize=10)
    a1[1,1].text(-1.5, -1.3, "Number of particles stagnating through the transition: " + str(trans), fontsize = 10)
    plt.savefig(f"{plot_loc}/stagnation_PT.png", dpi=500)
    print("plotted stagnant particles")

if __name__ == '__main__':
    main()
