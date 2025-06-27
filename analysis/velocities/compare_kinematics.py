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
pd.options.mode.chained_assignment = None  # default='warn'
import argparse
from pathlib import Path
import seaborn as sns
import matplotlib.tri as tri
path = str(Path(Path(__file__).parent.absolute()).parent.absolute())
sys.path.insert(0, path)
from libraries.functions import *
from libraries.particles import *
import plotly.graph_objects as go
from itertools import takewhile


def main():
    # Parse the arguments
    parser = argparse.ArgumentParser(description='Script that gets some models and gives the kinematic indicators')
    parser.add_argument('json_file', help='json file with model name, time at the end of computation, folder where to save the plots, legend')
    args = parser.parse_args()

    json_loc = '/home/vturino/PhD/projects/exhumation/pyInput/'
    json_file = f"{json_loc}{args.json_file}"

    # Read the json file
    with open(f"{json_file}") as json_file:
        configs = json.load(json_file)

    plotloc = f'/home/vturino/PhD/projects/exhumation/plots/{configs["plot_folder"][0]}'
    plotname = "kinematics_comparison.pdf"

    crfolder = f"/home/vturino/PhD/projects/exhumation/plots/single_models/{configs['models'][0]}/txt_files"

    cr = pd.read_csv(f"{crfolder}/2D_v.txt", sep="\s+")
    cr["conv_rate"].iloc[0] = np.nan

    ref_folder = f"/home/vturino/PhD/projects/exhumation/plots/single_models/kinematic_mu0.13_basalt7.5km_sed1km_cttV/txt_files"
    ref_dip = pd.read_csv(f"{ref_folder}/dip_data.csv")
    ref_trench = pd.read_csv(f"{ref_folder}/trench_pos.csv")


    # colors
    color = ["mediumslateblue", "darkslateblue", "mediumseagreen", "darkgreen"]
    markers = ["o", "s", "D", "v"]


    fig, ax = plt.subplots(1, 3, figsize=(15, 5))


    for m, mod in enumerate(configs['models']):
        mod_loc = f"/home/vturino/PhD/projects/exhumation/plots/single_models/{mod}"
        kin_loc = f"{mod_loc}/txt_files"
        if not os.path.exists(kin_loc):
            os.mkdir(kin_loc)
        dip = pd.read_csv(f"{kin_loc}/dip_data.csv")
        trench = pd.read_csv(f"{kin_loc}/trench_pos.csv")

        ax[0].plot(trench["time"], trench["trench_x"] - trench["trench_x"].iloc[0], color=color[m], label=configs['names'][m], linewidth=1.5)
        # ax[0].axhline(y=0, color='grey', linestyle='--', linewidth=1)
        ax[0].set_xlabel("Time (Ma)")
        ax[0].set_ylabel("Trench retreat (km)")
        ax[0].set_xlim(0, 50)


        ax[1].plot(dip["time"], dip["dip_shallow"], color=color[m], label=configs['names'][m], linewidth=1.5)  
        ax[1].set_xlabel("Time (Ma)")
        ax[1].set_ylabel("Shallow dip angle (degrees)")
        ax[1].set_xlim(0, 50)
        # ax[1].legend(loc='upper right', fontsize=10)

        ax[2].plot(dip["time"], dip["dip_deep"], color=color[m], label=configs['names'][m], linewidth=1.5)
        ax[2].set_xlabel("Time (Ma)")
        ax[2].set_ylabel("Deep dip angle (degrees)")
        ax[2].set_xlim(0, 50)
        # ax[2].legend(loc='upper right', fontsize=10)
       

    ax[0].plot(ref_trench["time"], ref_trench["trench_x"] - ref_trench["trench_x"].iloc[0], color="black", label="Reference model", linewidth=1.5, linestyle='--')
    ax[1].plot(ref_dip["time"], ref_dip["dip_shallow"], color="black", label="Reference model", linewidth=1.5, linestyle='--')
    ax[2].plot(ref_dip["time"], ref_dip["dip_deep"], color="black", label="Reference model", linewidth=1.5, linestyle='--')
    ax[0].legend(loc='upper right', fontsize=10)

    ax[0].set_title("Trench motion")
    ax[1].set_title("Shallow dip angle (50-100 km)")
    ax[2].set_title("Deep dip angle (150-300 km)")
    ax[0].grid(True)
    ax[1].grid(True)
    ax[2].grid(True)
    plt.tight_layout()
    
    plt.savefig(f"{plotloc}/{plotname}", dpi=300, bbox_inches='tight')
    plt.close("all")


    






if __name__ == "__main__":
    main()