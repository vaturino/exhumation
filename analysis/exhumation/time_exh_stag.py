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

    json_loc = '/home/vturino/PhD/projects/exhumation/pyInput/'
    # Read the json file
    with open(f"{json_loc}{args.json_file}") as json_file:
            configs = json.load(json_file)


    # Create the folders to save the plots
    plot_loc = f"/home/vturino/PhD/projects/exhumation/plots/single_models/{configs['models'][0]}"
    txt_loc = f'{plot_loc}/txt_files'
    if not os.path.exists(txt_loc):
        os.mkdir(txt_loc)

    # Load the data
    exhumed_list = pd.read_csv(f"{txt_loc}/timing_exhumed_particles.txt", sep="\s+")
    stagnant_list = pd.read_csv(f"{txt_loc}/timing_stagnant_particles.txt", sep="\s+")
    cr = pd.read_csv(f"{txt_loc}/2D_v.txt", sep="\s+")
    cr["conv_rate"].iloc[0]= np.nan

    #plot 2 rows of 3 plots, where we have, for each model, the exhumed and stagnant particles timing (tin, tmax, tfin) and the number of particles as a histogram.
    fig, axs = plt.subplots(2, 3, figsize=(15, 10))
    fig.suptitle(f"Timing of exhumed and stagnant particles for model {configs['models'][0]}")

    nbins = 10

    # Exhumed particles
    sns.histplot(data=exhumed_list, x="tin", bins=nbins, hue="lithology", palette="colorblind", alpha=1, ax=axs[0, 0], hue_order=exhumed_list["lithology"].value_counts(ascending=True).index)
    sns.histplot(data=exhumed_list, x="tmax", bins=nbins, hue="lithology", palette="colorblind", alpha=1, ax=axs[0, 1], hue_order=exhumed_list["lithology"].value_counts(ascending=True).index)   
    sns.histplot(data=exhumed_list, x="tfin", bins=nbins, hue="lithology", palette="colorblind", alpha=1, ax=axs[0, 2], hue_order=exhumed_list["lithology"].value_counts(ascending=True).index)
    sns.histplot(data=stagnant_list, x="tin", bins=nbins, hue="lithology", palette="colorblind", alpha=1, ax=axs[1, 0], hue_order=stagnant_list["lithology"].value_counts(ascending=True).index)
    sns.histplot(data=stagnant_list, x="tmax", bins=nbins, hue="lithology", palette="colorblind", alpha=1, ax=axs[1, 1], hue_order=stagnant_list["lithology"].value_counts(ascending=True).index)
    sns.histplot(data=stagnant_list, x="tfin", bins=nbins, hue="lithology", palette="colorblind", alpha=1, ax=axs[1, 2], hue_order=stagnant_list["lithology"].value_counts(ascending=True).index)
    

    
    ax1 = axs[0,0].twinx()
    ax1.plot(cr["time"]/1e6, cr["conv_rate"], color="grey", linewidth=1, zorder =10)
    ax1.set_ylabel("Convergence rate (cm/yr)", color="grey", fontweight='bold')
    axs[0,0].patch.set_visible(False) 
    axs[0,0].set_zorder(1) 
    ax1.set_zorder(10)
    for artist in axs[0, 0].get_children():
        if isinstance(artist, plt.Line2D):  # Ensure only lines are affected
            artist.set_zorder(0)  # Move histogram lines behind  
    
    ax2 = axs[0,1].twinx()
    ax2.plot(cr["time"]/1e6, cr["conv_rate"], color="grey", linewidth=2, zorder =10)
    ax2.set_ylabel("Convergence rate (cm/yr)", color="grey", fontweight='bold')
    axs[0,1].patch.set_visible(False)
    axs[0,1].set_zorder(1)
    ax2.set_zorder(10)
    for artist in axs[0, 1].get_children():
        if isinstance(artist, plt.Line2D):
            artist.set_zorder(0)

    ax3 = axs[0,2].twinx()
    ax3.plot(cr["time"]/1e6, cr["conv_rate"], color="grey", linewidth=2, zorder =10)
    ax3.set_ylabel("Convergence rate (cm/yr)", color="grey", fontweight='bold')
    axs[0,2].patch.set_visible(False)
    axs[0,2].set_zorder(1)
    ax3.set_zorder(10)
    for artist in axs[0, 2].get_children():
        if isinstance(artist, plt.Line2D):
            artist.set_zorder(0)

    ax4 = axs[1,0].twinx()
    ax4.plot(cr["time"]/1e6, cr["conv_rate"], color="grey", linewidth=2, zorder =10)
    ax4.set_ylabel("Convergence rate (cm/yr)", color="grey", fontweight='bold')
    axs[1,0].patch.set_visible(False)
    axs[1,0].set_zorder(1)
    ax4.set_zorder(10)
    for artist in axs[1, 0].get_children():
        if isinstance(artist, plt.Line2D):
            artist.set_zorder(0)

    ax5 = axs[1,1].twinx()
    ax5.plot(cr["time"]/1e6, cr["conv_rate"], color="grey", linewidth=2, zorder =10)
    ax5.set_ylabel("Convergence rate (cm/yr)", color="grey", fontweight='bold')
    axs[1,1].patch.set_visible(False)
    axs[1,1].set_zorder(1)
    ax5.set_zorder(10)
    for artist in axs[1, 1].get_children():
        if isinstance(artist, plt.Line2D):
            artist.set_zorder(0)

    ax6 = axs[1,2].twinx()
    ax6.plot(cr["time"]/1e6, cr["conv_rate"], color="grey", linewidth=2, zorder =10)
    ax6.set_ylabel("Convergence rate (cm/yr)", color="grey", fontweight='bold')
    axs[1,2].patch.set_visible(False)
    axs[1,2].set_zorder(1)
    ax6.set_zorder(10)
    for artist in axs[1, 2].get_children():
        if isinstance(artist, plt.Line2D):
            artist.set_zorder(0)


    plt.tight_layout()
    plt.savefig(f"{plot_loc}/timeline_exhumed_stagnant_particles.png", dpi = 500)
    plt.close()

if __name__ == "__main__":
    main()  