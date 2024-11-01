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
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.lines import Line2D


def main():
    parser = argparse.ArgumentParser(description= 'Script that gets some models and gives the kinematic indicators')
    parser.add_argument('json_file', help='json file with model name, time at the end of computation, folder where to save the plots, legend')
    args = parser.parse_args()

    json_loc = '/home/vturino/PhD/projects/exhumation/pyInput/'


    # Read the json file
    with open(f"{json_loc}{args.json_file}") as json_file:
            configs = json.load(json_file)
    m = configs["models"][0]

    # Create the folders to save the plots
    plot_loc = f"/home/vturino/PhD/projects/exhumation/plots/single_models/{m}"
    txt_loc = f'{plot_loc}/txt_files'
    if not os.path.exists(txt_loc):
        os.mkdir(txt_loc)


    # Load the data
    data = pd.read_csv(f"{txt_loc}/stagnant_maxima.txt", sep="\s+")
    cr = pd.read_csv(f"{txt_loc}/2D_v.txt", sep="\s+")
    cr["conv_rate"].iloc[0] = np.nan
    
    lnumber = len(data.lithology.unique())
    fig, ax = plt.subplots(2, lnumber, figsize=(17, 10)) 

    size = {
        "subducting": 100,
        "stagnant": 50,
        "rising": 20
    }

    color = {
        "Pmax": "#377eb8",  # Blue
        "P1": "#e41a1c",    # Red
        "P2": "#4daf4a",    # Green
        "P3": "#984ea3"     # Purple
    }

    colors_tmax = {
            "sed": "mediumblue",
            "oc": "#B06D1A",
            "ecl": "#45701C",
            "serp": "brown"
        }

    # Loop over each unique lithology to create separate subplots
    for ind_l, l in enumerate(data['lithology'].unique()):
        # Filter data and size values for each lithology
        lithology_data = data[data["lithology"] == l]

        
        # Plot Tmax vs Pmax for the current lithology
        sns.scatterplot(
            data=lithology_data, 
            x="Tmax", 
            y="Pmax", 
            ax=ax[0, ind_l], 
            size = "direction", sizes=size,
            legend = True,
            color = colors_tmax[l]
        )

    
        # # Customize each axis
        ax[0, ind_l].set_xlabel("Temperature (°C)", fontsize=13)
        ax[0, ind_l].set_ylabel("Peak Pressure (GPa)", fontsize=13)
        ax[0, ind_l].tick_params(axis='both', which='major', labelsize=12)
        ax[0, ind_l].set_ylim(0, data["Pmax"].max() + 0.05)
        ax[0, ind_l].set_xlim(0, data["Tmax"].max() + 50)
        ax[0, ind_l].set_title(l, fontsize=18, fontweight='bold', color = colors_tmax[l])
        
        # Plot time vs. pressure data in the second row of subplots
        sns.scatterplot(
            data=lithology_data, 
            x="timePmax", 
            y="Pmax", 
            size = "direction", sizes=size,
            ax=ax[1, ind_l],
            legend = False,
            color = colors_tmax[l],
        )
        ax1 = ax[1, ind_l].twinx()
        ax1.plot(cr["time"]/1e6, cr["conv_rate"], color="grey", linewidth=2, zorder = 10)
        ax1.set_ylabel("Convergence rate (cm/yr)", color="grey", fontweight='bold', fontsize=13)
        ax1.set_ylim(0, cr["conv_rate"].max() + 0.1)
        

            
        # Customize the time-pressure axes
        ax[1, ind_l].set_xlabel("Time (Myr)", fontsize=13)
        ax[1, ind_l].set_ylabel("Peak Pressure (GPa)", fontsize=13)
        ax[1, ind_l].tick_params(axis='both', which='major', labelsize=12)
        ax[1, ind_l].set_ylim(0, data["Pmax"].max() + 0.05)
        ax[1, ind_l].set_xlim(0, 50)


        fig.tight_layout()
        fig.savefig(f"{plot_loc}/stagnant_maxima.png", dpi=300)
        plt.close()

        # # Loop to plot additional (T1, P1), (T2, P2), (T3, P3) points for each lithology
        # for i in range(3):
        #     sns.scatterplot(
        #         data=lithology_data, 
        #         x=f"T{i+1}", 
        #         y=f"P{i+1}", 
        #         ax=ax[0, ind_l], 
        #         size = "direction", sizes=size,
        #         color = color[f"P{i+1}"],
        #         legend= False
        #     )

        # # Create a custom legend for size
        # size_legend = [Line2D([0], [0], marker='o', color='w', label=key, markersize=size[key] / 10, markerfacecolor='grey') 
        #             for key in size]
        # size_labels = list(size.keys())

        # # Create a custom legend for color
        # color_legend = [Line2D([0], [0], marker='o', color='w', label=key, markersize=10, markerfacecolor=color[key]) 
        #                 for key in color]
        # color_labels = list(color.keys())

        # # Place legends on the last subplot
        # ax[0, 0].legend(handles=size_legend, labels=size_labels, title="Direction", loc="upper left")
        # ax[0, 1].legend(handles=color_legend, labels=color_labels, title="Timing", loc="upper left")
        # # Loop to plot additional (timeP1, P1), (timeP2, P2), (timeP3, P3) points for each lithology
        # for i in range(3):
        #     sns.scatterplot(
        #         data=lithology_data, 
        #         x=f"timeP{i+1}", 
        #         y=f"P{i+1}", 
        #         ax=ax[1, ind_l],
        #         color = color[f"P{i+1}"],
        #         legend= False
        #     )


        




    






if __name__ == '__main__':
    main()