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

    compositions = configs['compositions']
    cutoff = configs['cutoff']

    file_count = 0

    for ind_m, m in tqdm(enumerate(configs['models'])):    
        time_array = np.zeros((len(os.listdir(f"{csvs_loc}{m}/fields")),2)) 
        stat = pd.read_csv(f"{models_loc}{m}/statistics",skiprows=configs['head_lines'],sep='\s+',header=None)
        time_array = grab_dimTime_fields(f"{csvs_loc}{m}/fields", stat, time_array, configs['head_lines']-1)
        plot_loc_mod = f"/home/vturino/PhD/projects/exhumation/plots/single_models/{m}"
        if not os.path.exists(plot_loc_mod):
            os.mkdir(plot_loc_mod)
        plot_loc = f"{plot_loc_mod}/Panels/"
        if not os.path.exists(plot_loc):
            os.mkdir(plot_loc)
        else:
            # Count the files in the fields_loc directory
            file_count = len(os.listdir(plot_loc))


        for t in tqdm(range(2*file_count, len(time_array), 2)):

            f1, a1 =plt.subplots(2, 2, figsize=(15, 12))
            plotname = f"{plot_loc}{int(t/2)}.png" 
            data = pd.read_parquet(f"{csvs_loc}{m}/fields/full.{int(t)}.gzip") 
            data["lithology"] = 0
            data["logvisc"] = np.log10(data["viscosity"])
            data["logSR"] = np.log10(data["strain_rate"])

            # create lithology classes
            for ind, c in enumerate(compositions):
                data[c][data["Points:1"] < cutoff[ind]] = 0
                data[c][data[c] >= 0.5] = 1
                data[c][data[c] < 0.5] = 0

            for ind_c, c in enumerate(compositions):
                weight = ind_c + 1 
                data["lithology"] += weight * data[c]
                
            data["lithology"] = data["lithology"].astype(int)
            
            
            
            
            pts = get_points_with_y_in(data, 15.e3, 2.e3, ymax = 900.e3)
            trench= get_trench_position(pts,threshold = 0.13e7)
            xmin_plot = trench -100.e3
            xmax_plot = trench + 200.e3
            ymin_plot = 740.e3
            ymax_plot = 902.e3

            x = data["Points:0"].to_numpy()/1.e3
            y = (ymax_plot - data["Points:1"].to_numpy())/1.e3

            triang = tri.Triangulation(x, y)
            colors = matplotlib.colormaps['Accent'].resampled(len(compositions)+1)
            colors.set_bad('white')
           
            # plot only triangles with sidelength smaller some max_radius
            max_radius = 10
            triangles = triang.triangles

            # Mask off unwanted triangles.  
            xtri = x[triangles] - np.roll(x[triangles], 1, axis=1)
            ytri = y[triangles] - np.roll(y[triangles], 1, axis=1)
            maxi = np.max(np.sqrt(xtri**2 + ytri**2), axis=1)
            triang.set_mask(np.logical_and(maxi > max_radius, y[triangles][:,1] < 90.))

            
            # row = exhumed_data[exhumed_data['time'] == t]
            
            
            p1 = a1[0,0].tripcolor(triang, data["lithology"], cmap=colors, shading='gouraud', vmax=len(compositions), vmin=0)
            a1[0,0].spines[['top']].set_visible(False)
            plt.colorbar(p1, orientation='horizontal', label='Composition', ax = a1[0,0])
            a1[0,0].set_ylim([(ymax_plot-ymin_plot)/1.e3,-5])
            a1[0,0].set_xlim([xmin_plot/1.e3,xmax_plot/1.e3])
            a1[0,0].set_aspect('equal')

            p2 = a1[0,1].tripcolor(triang, data["logSR"], cmap='RdBu_r', shading='gouraud', vmin=-19, vmax=-12)
            a1[0,1].spines[['top']].set_visible(False)
            plt.colorbar(p2, orientation='horizontal', label='Log(Strain rate) [s-1]', ax = a1[0,1])
            a1[0,1].set_ylim([(ymax_plot-ymin_plot)/1.e3,-5])
            a1[0,1].set_xlim([xmin_plot/1.e3,xmax_plot/1.e3])
            a1[0,1].set_aspect('equal')

            p3 = a1[1,0].tripcolor(triang, data["logvisc"], shading='gouraud', vmin=18, vmax=24)
            a1[1,0].spines[['top']].set_visible(False)
            plt.colorbar(p3, orientation='horizontal', label='Log(Viscosity) [Pa s]', ax = a1[1,0])
            a1[1,0].set_ylim([(ymax_plot-ymin_plot)/1.e3,-5])
            a1[1,0].set_xlim([xmin_plot/1.e3,xmax_plot/1.e3])
            a1[1,0].set_aspect('equal')

            p4 = a1[1,1].tripcolor(triang, data["T"]-273.5, cmap='inferno', shading='gouraud', vmin=0, vmax=1421)
            a1[1,1].tricontour(triang, data["T"]-273.5, colors='white', levels=[100, 300, 500, 700, 900], linewidths=0.5)
            a1[1,1].spines[['top']].set_visible(False)
            plt.colorbar(p4, orientation='horizontal', label='Temperature [C]', ax = a1[1,1])
            a1[1,1].set_ylim([(ymax_plot-ymin_plot)/1.e3,-5])
            a1[1,1].set_xlim([xmin_plot/1.e3,xmax_plot/1.e3])
            a1[1,1].set_aspect('equal')

            plt.savefig(plotname, bbox_inches='tight', format='png', dpi=1000)
            plt.clf()
            plt.close('all')

if __name__ == "__main__":
    main()

