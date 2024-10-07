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


## FUNCTIONS
def grab_dimTime_fields (dir: str, stat:str, time, skiprows = 15):
    ts = len(os.listdir(dir))
    for t in range(ts):
        time[t,0] = t
    filt = stat.loc[:,skiprows].notnull()
    time[:,1] = stat[filt].loc[:,1]
    return time




def main():
    parser = argparse.ArgumentParser(description= 'Script that gets n models and the time at the end of computation and gives the temperature and viscosity plots for all time steps')
    parser.add_argument('json_file', help='json file with models and end time to plot T and eta field')  
    args = parser.parse_args()

    json_loc = 'path to jason file with model name, plot folder and header lines for statistics file'

    models_loc = 'where the raw output is'
    csvs_loc = 'where the plots are'

    with open(f"{json_loc}{args.json_file}") as json_file:
            configs = json.load(json_file)


    xmin_plot = 0
    xmax_plot = 5400.e3
    ymin_plot = 0
    ymax_plot = 905.e3


    for ind_m, m in tqdm(enumerate(configs['models'])):    
        time_array = np.zeros((len(os.listdir(f"{csvs_loc}{m}/fields")),2)) 
        stat = pd.read_csv(f"{models_loc}{m}/statistics",skiprows=configs['head_lines'],sep='\s+',header=None)
        time_array = grab_dimTime_fields(f"{csvs_loc}{m}/fields", stat, time_array, configs['head_lines']-1)
        plot_loc_mod = f"/home/vturino/PhD/projects/exhumation/plots/single_models/{m}"
        if not os.path.exists(plot_loc_mod):
            os.mkdir(plot_loc_mod)
        plot_loc = f"{plot_loc_mod}/Viscosity/"
        if not os.path.exists(plot_loc):
            os.mkdir(plot_loc)

        
        for t in tqdm(range(0, len(time_array), 2)):
            plotname = f"{plot_loc}{int(t/2)}.png" 
            data = pd.read_parquet(f"{csvs_loc}{m}/fields/full.{int(t)}.gzip") 
            data["logvisc"] = np.log10(data["viscosity"])

            # Convert x, y from data into arrays so we can build a mesh
            x = data["Points:0"].to_numpy()/1.e3
            y = (ymax_plot - data["Points:1"].to_numpy())/1.e3

            # Create a meshgrid with triangulation function
            triang = tri.Triangulation(x, y)
           
           
            # plot only triangles with sidelength smaller some max_radius
            max_radius = 10
            triangles = triang.triangles

            # Mask off unwanted triangles.  
            xtri = x[triangles] - np.roll(x[triangles], 1, axis=1)
            ytri = y[triangles] - np.roll(y[triangles], 1, axis=1)
            maxi = np.max(np.sqrt(xtri**2 + ytri**2), axis=1)
            triang.set_mask(np.logical_and(maxi > max_radius, y[triangles][:,1] < 90.))

            # Plot
            plt.tripcolor(triang, data["logvisc"], shading='gouraud', vmin=18, vmax=24)
            plt.colorbar(orientation='horizontal', label='Log(Viscosity) [Pa s]')
            plt.ylim([(ymax_plot-ymin_plot)/1.e3,-5])
            plt.xlim([xmin_plot/1.e3,xmax_plot/1.e3])
            plt.gca().set_aspect('equal')

            plt.savefig(plotname, bbox_inches='tight', format='png', dpi=1000)
            plt.clf()
            plt.close('all')

if __name__ == "__main__":
    main()

