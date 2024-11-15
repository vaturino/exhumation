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

    # Load the data
    allp = pd.read_csv(f"{txt_loc}/particles_indexes.txt", sep="\s+")
    stag = pd.read_csv(f"{txt_loc}/stagnant_particles.txt", sep="\s+")
    exh = pd.read_csv(f"{txt_loc}/exhumed_particles.txt", sep="\s+")
    subd = pd.read_csv(f"{txt_loc}/subducted_particles.txt", sep="\s+")
    subd = subd.dropna()

    # Check unique values in stag['lithology']
    unique_lithologies_stag = stag["lithology"].unique()
    unique_lithologies_exh = exh["lithology"].unique()

    # Create a dictionary mapping each unique lithology to its corresponding color
    palette_stag = {lith: colors_tfin[lith] for lith in unique_lithologies_stag if lith in colors_tfin}
    palette_exh = {lith: colors_tfin[lith] for lith in unique_lithologies_exh if lith in colors_tfin}

    # Sample counts for subducted, exhumed, and stagnant particles
    subd_count = len(subd)
    exh_count = len(exh)
    stag_count = len(stag)

    # Total count of all particles
    total_count = subd_count + exh_count + stag_count

    # Calculate percentages for each category
    subd_pct = (subd_count / total_count) * 100
    exh_pct = (exh_count / total_count) * 100
    stag_pct = (stag_count / total_count) * 100


    #calculate percentages of exhumed particles by lithology
    exh_pct_lith = exh.groupby("lithology").size() / total_count * 100
    stag_pct_lith = stag.groupby("lithology").size() / total_count * 100


    bar_fontsize = 11


    


    f1, a1 = plt.subplots(1, 4, figsize=(17, 5))

    #plot exhumed peak P conditions over rocks
    sns.scatterplot(data=exh, 
                    x="maxPT", 
                    y="maxPP", 
                    hue="lithology",  
                    ax=a1[1], 
                    zorder = 10,
                    alpha=1,
                    palette=palette_exh)
    a1[1].set_ylabel("Pressure (GPa)")
    a1[1].set_xlabel("T ($^\circ$C)")
    a1[1].set_title("Exhumed particles: peak Pressure")
    a1[1].set_xlim(0, 900)
    a1[1].set_ylim(0, 3.5)
    sns.kdeplot(data=rocks, x="T", y="P", ax=a1[1], fill = True, cbar = False, alpha = .7, zorder=2, color='grey')
    

    #plot stagnating peak P conditions over rocks
    sns.scatterplot(data=stag, x="Tm_kin", y="Pm_kin", hue="lithology", ax=a1[2], zorder = 10, palette=palette_stag, legend=True)
    sns.scatterplot(data=stag, x="Tm_dyn", y="Pm_dyn", hue="lithology", ax=a1[2], zorder = 10, palette=palette_stag, legend=False)
    sns.scatterplot(data=stag, x="Tm_trans", y = "Pm_trans", hue="lithology", ax=a1[2], zorder = 10, palette=palette_stag, legend=False)
    a1[2].set_ylabel("Pressure (GPa)")
    a1[2].set_xlabel("T ($^\circ$C)")
    a1[2].set_title("Peak pressure")
    a1[2].set_xlim(0, 900)
    a1[2].set_ylim(0, 3.5)
    sns.kdeplot(data=rocks, x="T", y="P", ax=a1[2], fill = True, cbar = False, alpha = .7, zorder=2, color='grey')
    
    


    for s in range(0, len(subd), 50):
        id = subd["id"].iloc[s].astype(int)
        spart = pd.read_csv(f"{pt_files}/pt_part_{id}.txt", sep="\s+")
        spart["Plith"] = (900.- spart["depth"])*1e3*9.81*3100/1e9
        a1[0].plot((spart["T"]), spart["Plith"], color = 'darkslateblue', alpha = 0.5, linewidth = 0.5, zorder = 20)
    a1[0].set_ylabel("Pressure (GPa)")
    a1[0].set_xlabel("Temperature (C)")
    a1[0].set_title("Subducted particles")
    sns.kdeplot(data=rocks, x="T", y="P", ax=a1[0], fill = True, cbar = False, alpha = .7, zorder=1, color='grey')
    a1[0].set_ylim(0, 3.5)
    a1[0].set_xlim(0, 900)



   # Plot the 100% stacked vertical bar with a narrow width for Particle Classification
    a1[3].bar([0], subd_pct, color='powderblue', label='Subducted', width=0.1, linewidth=0.5)
    a1[3].bar([0], exh_pct, bottom=subd_pct, color='cornflowerblue', label='Exhumed', width=0.1)
    a1[3].bar([0], stag_pct, bottom=subd_pct + exh_pct, color='olivedrab', label='Stagnant', width=0.1)

    # Set limits and formatting for the particle classification bar
    a1[3].set_ylim(0, 100)
    a1[3].set_xlim(-0.15, 0.85)  # Adjust x-limits to accommodate all the bars
    a1[3].set_xticks([0])  # Set x-ticks to only display the "Total" position
    a1[3].set_xticklabels(['Total'])  # Label the x-tick as 'Total'
    a1[3].set_yticks([])
    a1[3].set_title('Particle Classification')

    # Add percentages inside the bars
    a1[3].text(0, subd_pct / 2, f'{subd_pct:.1f}%', ha='center', va='center', color='k', fontsize=bar_fontsize)
    a1[3].text(0, subd_pct + exh_pct / 2, f'{exh_pct:.1f}%', ha='center', va='center', color='k', fontsize=bar_fontsize)
    a1[3].text(0, subd_pct + exh_pct + stag_pct / 2, f'{stag_pct:.1f}%', ha='center', va='center', color='k', fontsize=bar_fontsize)
    a1[3].spines[['top', 'right', 'left', 'bottom']].set_visible(False)
    a1[3].legend(loc='upper left')

    # Plot the 100% stacked vertical bar for lithology distribution in Exhumed particles
    bottom = 0
    for lith, pct in exh_pct_lith.items():
        a1[3].bar([0.2], pct * 4, bottom=bottom, color=palette_exh[lith], label=lith, width=0.1)
        a1[3].text(0.2, bottom + (pct * 2), f'{pct:.1f}%', ha='center', va='center', color='k', fontsize=bar_fontsize)
        bottom += pct * 4

    # Plot the 100% stacked vertical bar for lithology distribution in Stagnant particles
    bottom = 0
    for lith, pct in stag_pct_lith.items():
        a1[3].bar([0.4], pct * 4, bottom=bottom, color=palette_stag[lith], width=0.1)
        a1[3].text(0.4, bottom + (pct * 2), f'{pct:.1f}%', ha='center', va='center', color='k', fontsize=bar_fontsize)
        bottom += pct * 4

    # Set limits and formatting for lithology distribution bars
    a1[3].set_ylim(0, 100)
    a1[3].set_xlim(-0.15, 0.65)  # Adjust x-limits to accommodate all the bars
    a1[3].set_xticks([0, 0.2, 0.4])  # Set x-ticks to center each lithology bar
    a1[3].set_xticklabels(['Total', 'Exhumed', 'Stagnant'])  # Label each group
    a1[3].set_yticks([])  # Optional: fewer y-ticks to declutter

    # Add the legend for lithology distribution
    a1[3].legend(loc='upper center')

    

    

        
    



    

    f1.tight_layout()
    
    format = ["png", "eps"]
    for f in format:
        plt.savefig(f"{plot_loc}/peak_P_comparison.{f}", dpi = 500)
    # plt.savefig(f"{plot_loc}/percentages.png", dpi = 500)
    plt.close()
         



if __name__ == '__main__':
    main()