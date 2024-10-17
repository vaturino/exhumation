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
        plot_loc = f"{plot_loc_mod}/Zoomed_visc/"
        if not os.path.exists(plot_loc):
            os.mkdir(plot_loc)
        else:
            # Count the files in the fields_loc directory
            file_count = len(os.listdir(plot_loc))

        compositions = configs['compositions']
        cutoff = configs['cutoff']


        
        for t in tqdm(range(0, len(time_array), 2)):
        # for t in tqdm(range(2,3)):

            fig=plt.figure()
            gs=GridSpec(2,1)
            plotname = f"{plot_loc}{int(t/2)}.png" 
            data = pd.read_parquet(f"{csvs_loc}{m}/fields/full.{int(t)}.gzip") 
            data["lithology"] = 0
            data["logvisc"] = np.log10(data["viscosity"])
            for ind, c in enumerate(compositions):
                data[c][data["Points:1"] < cutoff[ind]] = 0
                data[c][data[c] >= 0.5] = 1
                data[c][data[c] < 0.5] = 0

            for ind_c, c in enumerate(compositions):
                weight = ind_c + 1
                data["lithology"] += weight * data[c]
            
            data["lithology"] = data["lithology"].astype(int)
            composition_mapping = {ind_c + 1: c for ind_c, c in enumerate(compositions)}
            data["terrain"] = data["lithology"].map(composition_mapping)
            data["terrain"].fillna("mantle", inplace=True)
            filter_mask = (data["terrain"] != "mantle") & (data["terrain"] != "opc")
            data_filtered_comp = data[filter_mask]

            
            pts = get_points_with_y_in(data, 15.e3, 2.e3, ymax = 900.e3)
            trench= get_trench_position(pts,threshold = 0.13e7)
            xmin_plot = trench -100.e3
            xmax_plot = trench + 200.e3
            ymin_plot = 740.e3
            ymax_plot = 900.e3

            x = data["Points:0"].to_numpy()/1.e3
            y = (ymax_plot - data["Points:1"].to_numpy())/1.e3

            triang = tri.Triangulation(x, y)
            colors = matplotlib.colormaps['viridis']((18,24,601))
           
           
            # plot only triangles with sidelength smaller some max_radius
            max_radius = 10
            triangles = triang.triangles

            # Mask off unwanted triangles.  
            xtri = x[triangles] - np.roll(x[triangles], 1, axis=1)
            ytri = y[triangles] - np.roll(y[triangles], 1, axis=1)
            maxi = np.max(np.sqrt(xtri**2 + ytri**2), axis=1)
            triang.set_mask(np.logical_and(maxi > max_radius, y[triangles][:,1] < 90.))

          
            plt.tripcolor(triang, data["logvisc"], shading='gouraud', vmin=18, vmax=24)
            vel_plot_thresh = 0.0001  # don't plot velocity vectors smaller than this (cm/yr)
            step = 500  # plot every 100th vector to reduce the number of vectors plotted
            # Filter out small velocity vectors and downsample data for plotting
            mask = (100. * np.sqrt(data["velocity:0"]**2 + data["velocity:1"]**2)) >= vel_plot_thresh
            mask = mask.reindex(data_filtered_comp.index, fill_value=False)
            data_filtered = data_filtered_comp[mask].iloc[::step]

            vel_vects = plt.quiver(data_filtered["Points:0"].to_numpy() / 1.e3, 
                                     (ymax_plot - data_filtered["Points:1"].to_numpy()) / 1.e3, 
                                     data_filtered["velocity:0"].to_numpy() * 100, 
                                     data_filtered["velocity:1"].to_numpy() * 100, 
                                     scale=50, color='black', width=0.0015)
            plt.quiverkey(vel_vects, 0.15, 0.1, 1, '1 cm/yr', labelpos='W', fontproperties={'size': '7'}, color='white', labelcolor='white')
            # plt.scatter(row["x"].to_numpy()/1.e3, (ymax_plot - row["depth"].to_numpy())/1.e3, color='red', marker='x', zorder = 1)
            plt.colorbar(orientation='horizontal', label='Log(Viscosity) [Pa s]')
            plt.ylim([(ymax_plot-ymin_plot)/1.e3,-5])
            plt.xlim([xmin_plot/1.e3,xmax_plot/1.e3])
            plt.gca().set_aspect('equal')
            # plt.xticks(np.arange(xmin_plot/1.e3, xmax_plot/1.e3, 20))
            # plt.yticks(np.arange((ymax_plot-ymin_plot)/1.e3, -5, -20))

            plt.savefig(plotname, bbox_inches='tight', format='png', dpi=1000)
            plt.clf()
            plt.close('all')

if __name__ == "__main__":
    main()

