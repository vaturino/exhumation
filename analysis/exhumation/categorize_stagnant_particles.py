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
from concurrent.futures import ProcessPoolExecutor, as_completed

path = str(Path(Path(__file__).parent.absolute()).parent.absolute())
sys.path.insert(0, path)
from libraries.functions import *
from libraries.particles import *
from libraries.exhumation import *  
import plotly.express as px





############### MAIN ####################

def main():
    parser = argparse.ArgumentParser(description='Script that gets some models and gives the kinematic indicators')
    parser.add_argument('json_file', help='json file with model name, time at the end of computation, folder where to save the plots, legend')
    args = parser.parse_args()

    json_loc = '/home/vturino/PhD/projects/exhumation/pyInput/'

    with open(f"{json_loc}{args.json_file}") as json_file:
        configs = json.load(json_file)
    m = configs["models"][0]

    plot_loc = f"/home/vturino/PhD/projects/exhumation/plots/single_models/{m}"
    txt_loc = f'{plot_loc}/txt_files'
    if not os.path.exists(txt_loc):
        os.mkdir(txt_loc)

    txt_loc = f'{plot_loc}/txt_files'

    # Read the data from the CSV file
    data = pd.read_csv(f"{txt_loc}/stagnant_particles.txt", sep="\s+")
    cr = pd.read_csv(f"{txt_loc}/2D_v.txt", sep="\s+")
    cr.iloc[0] = np.nan
    cr = cr.dropna()

    # one = data[data[['Pm_kin', 'Pm_dyn', 'Pm_trans']].notna().sum(axis=1) == 1]
    # multiple = data[~data.index.isin(data.index)]


    # f1, a1 = plt.subplots(2, 1, figsize=(15, 7), height_ratios=[0.25, 1])
    # size = 5
    # alpha_line = 0.5

    # # Plot convergence rate
    # a1[0].plot(cr['time']/1e6, cr['conv_rate'], color='grey', label='Convergence rate', linewidth=2)
    # a1[0].label_outer()  # Only show outer labels and tick labels
    # a1[0].set_xlim(0, 50)
    # a1[0].set_ylim(0, max(cr['conv_rate']) + 0.2)
    # a1[0].axvline(x=35., color='grey', linewidth=1, linestyle="--")
    
    

    # # Plot horizontal lines for each data point with range from ti_kin to ti_kin + time_interval_kin
    # for i in range(len(data)):
    #     ti_kin_value = data['ti_kin'].iloc[i]
    #     time_interval_kin_value = data['time_interval_kin'].iloc[i]
    #     Pm_kin_value = data['Pm_kin'].iloc[i]
        
    #     # Plot horizontal line from (tm_kin - time_interval_kin) to (tm_kin + time_interval_kin)
    #     a1[1].plot([ti_kin_value, ti_kin_value + time_interval_kin_value], 
    #             [Pm_kin_value, Pm_kin_value], color='lightblue', linewidth=0.5, alpha= alpha_line)
    # a1[1].scatter(data['tm_kin'], data['Pm_kin'], color='royalblue', label='Data points', zorder = 10, alpha = 1, s = size)

    # # Plot horizontal lines for each data point with range from tmi_dyn to tmi_dyn + time_interval_dyn
    # for i in range(len(data)):
    #     ti_dyn_value = data['ti_dyn'].iloc[i]
    #     time_interval_dyn_value = data['time_interval_dyn'].iloc[i]
    #     Pm_dyn_value = data['Pm_dyn'].iloc[i]
        
    #     # Plot horizontal line from (tm_dyn - time_interval_dyn) to (tm_dyn + time_interval_dyn)
    #     a1[1].plot([ti_dyn_value, ti_dyn_value + time_interval_dyn_value], 
    #             [Pm_dyn_value, Pm_dyn_value], color='lightcoral', linewidth=0.5, alpha= alpha_line)
    # a1[1].scatter(data['tm_dyn'], data['Pm_dyn'], color='firebrick', label='Data points', zorder = 10, alpha = 1, s = size)

    # # Plot horizontal lines for each data point with range from ti_trans to ti_trans+time_interval_trans
    # for i in range(len(data)):
    #     ti_trans_value = data['ti_trans'].iloc[i]
    #     time_interval_trans_value = data['time_interval_trans'].iloc[i]
    #     Pm_trans_value = data['Pm_trans'].iloc[i]
        
    #     # Plot horizontal line from (tm_dyn - time_interval_dyn) to (tm_dyn + time_interval_dyn)
    #     a1[1].plot([ti_trans_value, ti_trans_value + time_interval_trans_value], 
    #             [Pm_trans_value, Pm_trans_value], color='blueviolet', linewidth=1, alpha= alpha_line)
    # a1[1].scatter(data['tm_trans'], data['Pm_trans'], color='indigo', label='Data points', zorder = 10, alpha = 1, s = size)

    # # Add labels and title to the plot
    # a1[1].set_xlabel('time (Myr)')
    # a1[1].set_ylabel('Pressure')
    # a1[1].set_xlim(0, 50)
    # a1[1].axvline(x=35., color='grey', linewidth=1, linestyle="--")
    # # a1[1].set_title('Stagnation time range')
    # # a1[1].legend()

    

    # plt.subplots_adjust(hspace=0.)


    # plt.savefig(f"{plot_loc}/stagnant_time_intervals.png")
    # plt.close()



    #now I want to group by pressure: I want to have a list of particles that stagnate at a certain pressure interval, 
    # compute the average time of stagnation and the number of particles that stagnate at that pressure interval and 
    # the average intial time and plot them. Max 3 values of pressure for each group.
    # Do for each Pm_kin, Pm_dyn, Pm_trans

    # Group by pressure
    min = data['Pm_kin'].min()
    max = data['Pm_kin'].max()
    #Find 3 values in this interval
    step = (max - min)/3
    pressure_intervals = []
    for i in range(3):
        pressure_intervals.append((min + step*i, min + step*(i+1)))
    print(pressure_intervals)



   


    # plt.show()


    


if __name__ == "__main__":
    main()
