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
    parser = argparse.ArgumentParser(description= 'Script that gets n models and the time at the end of computation and gives the temperature and viscosity plots for all time steps')
    parser.add_argument('json_file', help='json file with models and end time to plot T and eta field')  
    args = parser.parse_args()

    # csvs_loc =  '/home/vturino/PhD/projects/exhumation/gz_outputs/'
    # models_loc =  '/home/vturino/PhD/projects/exhumation/raw_outputs/'
    json_loc = '/home/vturino/PhD/projects/exhumation/pyInput/'

    models_loc = '/home/vturino/Vale_nas/exhumation/raw_outputs/'
    csvs_loc = '/home/vturino/Vale_nas/exhumation/gz_outputs/'

    with open(f"{json_loc}{args.json_file}") as json_file:
            configs = json.load(json_file)

    # compositions = configs['compositions']
    # cutoff = configs['cutoff']

    file_count = 0


    for ind_m, m in tqdm(enumerate(configs['models'])):    
        time_array = np.zeros((len(os.listdir(f"{csvs_loc}{m}/fields")),2)) 
        stat = pd.read_csv(f"{models_loc}{m}/statistics",skiprows=configs['head_lines'],sep='\s+',header=None)
        time_array = grab_dimTime_fields(f"{csvs_loc}{m}/fields", stat, time_array, configs['head_lines']-1)
        plot_loc_mod = f"/home/vturino/PhD/projects/exhumation/plots/single_models/{m}"
        if not os.path.exists(plot_loc_mod):
            os.mkdir(plot_loc_mod)
        plot_loc = f"{plot_loc_mod}/kinematics/"
        if not os.path.exists(plot_loc):
            os.mkdir(plot_loc)
        else:
            # Count the files in the fields_loc directory
            file_count = len(os.listdir(plot_loc))

        points_loc = f"{plot_loc_mod}/dip_calculation_points/"
        if not os.path.exists(points_loc):
            os.mkdir(points_loc)
        else:
            # Count the files in the fields_loc directory
            file_count = len(os.listdir(points_loc))
        
        txt_loc = f"{plot_loc_mod}/txt_files/"
        if not os.path.exists(txt_loc):
            os.mkdir(txt_loc)
        else:
            # Count the files in the fields_loc directory
            file_count = len(os.listdir(txt_loc))

        
        xmin_plot = 0
        xmax_plot = 5400.e3
        ymin_plot = 0
        ymax_plot = 902.e3

        
        plotname = "dip.png"
        dip_data = pd.DataFrame(columns=["time", "dip_shallow", "x1", "y1", "x2", "y2", "dip_deep", "x3", "y3", "x4", "y4"], index=range(len(time_array)))
        
        for t in tqdm(range(0, len(time_array))):
        # for t in tqdm(range(0,5)):
            fig, ax = plt.subplots()
            data = pd.read_parquet(f"{csvs_loc}{m}/fields/full.{int(t)}.gzip") 
            data["dip"] = 0
            data["logvisc"] = np.log10(data["viscosity"])

            x = data["Points:0"].to_numpy()/1.e3
            y = (ymax_plot - data["Points:1"].to_numpy())/1.e3

            triang = tri.Triangulation(x, y)
    

            core_cont = ax.tricontour(triang, data["core"], levels=[0.5], colors='black', linewidths=0.5, linestyles='solid', alpha=0)
            core_cont_array = core_cont.allsegs[0][0]

            
            # just use top contour 
            max_y = np.max(core_cont_array[:, 1])
            max_y_index = np.where(core_cont_array[:, 1] == max_y)[0][0]
            top_contour = core_cont_array[:max_y_index + 1]
            # smooth out the contour
            top_contour[:, 0] = pd.Series(top_contour[:, 0]).rolling(window=10, center=True).mean().to_numpy()
            top_contour[:, 1] = pd.Series(top_contour[:, 1]).rolling(window=10, center=True).mean().to_numpy()
            mask = ~np.isnan(top_contour).any(axis=1)
            top_contour = top_contour[mask]


            flat_mask = top_contour[:, 1] <35.
            x_max_flat = np.max(top_contour[flat_mask, 0])
            index_flat = np.where(top_contour[:, 0] == x_max_flat)[0][0]
            
            x1, y1 = top_contour[index_flat, 0], top_contour[index_flat, 1]
            #y2 is y that is closest to 200
            target_y = 100
            diffs = np.abs(top_contour[:, 1] - target_y)
            index_y2 = np.argmin(diffs)
            x2, y2 = top_contour[index_y2, 0], top_contour[index_y2, 1]

            target_d1 = 150
            target_d2 = 300
            d1 = np.abs(top_contour[:,1] - target_d1)
            d2 = np.abs(top_contour[:, 1] - target_d2)
            index_d1 = np.argmin(d1)
            index_d2 = np.argmin(d2)
            x3, y3 = top_contour[index_d1, 0], top_contour[index_d1, 1]
            x4, y4 = top_contour[index_d2, 0], top_contour[index_d2, 1]

            # calculate dip angle
            dip_angle_shallow = np.arctan2(y2 - y1, x2 - x1) * 180 / np.pi

            
            dip_angle_deep = np.arctan2(y4 - y3, x4 - x3) * 180 / np.pi
            if dip_angle_deep > 90 or dip_angle_deep == 0:
                dip_angle_deep = np.nan
            dip_data.iloc[t] = [time_array[t, 1]/1.e6, dip_angle_shallow, x1, y1, x2, y2, dip_angle_deep, x3, y3, x4, y4]

            ax.plot(top_contour[:, 0], top_contour[:, 1], color='black', linewidth=0.5, label='Top Contour')
            ax.axhline(y=y1, color='red', linestyle='--', linewidth=0.5, alpha = 0.5)
            ax.plot([x1, x2], [y1, y2], color='red', linestyle='-', linewidth=0.5, alpha = 0.8, marker = "o", markersize=1.5, label='Shallow Dip')
            ax.axhline(y=y3, color='blue', linestyle='--', linewidth=0.5, alpha = 0.5)
            ax.plot([x3, x4], [y3, y4], color='blue', linestyle='-', linewidth=0.5, alpha = 0.8, marker = "o", markersize=1.5, label='Deep Dip')
            ax.set_ylim([(ymax_plot-ymin_plot)/1.e3,-5])
            ax.set_xlim([xmin_plot/1.e3,xmax_plot/1.e3])
            ax.spines['top'].set_visible(False)
            ax.set_aspect('equal')
            ax.set_xlabel("X (km)")
            ax.set_ylabel("Y (km)")
            ax.legend(loc='lower right')
            ax.set_title(f"Dip Angles at t={(time_array[t, 1]/1.e6):.1f} Ma: Shallow = {dip_angle_shallow:.2f}$^\circ$; Deep = {dip_angle_deep:.2f}$^\circ$")

            fig.tight_layout()
            fig.savefig(f"{points_loc}/dip_points_{int(t)}.png", bbox_inches='tight', format='png', dpi=500)

            plt.close('all')

    #save dip_data to csv
    dip_data.to_csv(f"{txt_loc}/dip_data.csv", index=False)

    plt.plot(dip_data["time"], dip_data["dip_shallow"], marker='o', label='Shallow Dip', color='red', markersize=2, linewidth=0.5, alpha=0.8)
    plt.plot(dip_data["time"], dip_data["dip_deep"], marker='o', label='Deep Dip', color='blue', markersize=2, linewidth=0.5, alpha=0.8)
    plt.legend()
    plt.xlabel("Time (s)")
    plt.ylabel("Dip Angle (degrees)")
    plt.title("Dip Angle Over Time")
    plt.savefig(f"{plot_loc}{plotname}", bbox_inches='tight', format='png', dpi=500)

    plt.close("all")






if __name__ == "__main__":
    main()