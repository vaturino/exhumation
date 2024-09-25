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

    # csvs_loc =  '/home/vturino/PhD/projects/exhumation/gz_outputs/'
    # models_loc =  '/home/vturino/PhD/projects/exhumation/raw_outputs/'

    models_loc = '/home/vturino/Vale_nas/exhumation/raw_outputs/'
    csvs_loc = '/home/vturino/Vale_nas/exhumation/gz_outputs/'

    json_loc = '/home/vturino/PhD/projects/exhumation/pyInput/'

    freq = 1


    with open(f"{json_loc}{args.json_file}") as json_file:
            configs = json.load(json_file)
            setting = 'none'

    compositions = configs['compositions']
    
    for ind_m, m in tqdm(enumerate(configs['models'])): 
        time_array = np.zeros((len(os.listdir(f"{csvs_loc}{m}/particles")),2))   
        stat = pd.read_csv(f"{models_loc}{m}/statistics",skiprows=configs['head_lines'],sep='\s+',header=None)
        time_array = grab_dimTime_particles(f"{csvs_loc}{m}/particles", stat, time_array, configs['head_lines']-1)
        ts = int(len(time_array))

        # trench = np.zeros(ts-1)

        # for t in range(1, ts):
        #     trench[t-1] = get_trench_position(pd.read_parquet(f"{csvs_loc}{m}/fields/full.{t-1}.gzip"), threshold = 1.3e6)


        plot_loc = f"/home/vturino/PhD/projects/exhumation/plots/single_models/{configs['models'][ind_m]}"
        # plot_loc_init = f"/home/vturino/PhD/projects/exhumation/plots/single_models/kinematic_mu0.025_upto35Myr"
        txt_loc = f'{plot_loc}/txt_files'
        if not os.path.exists(txt_loc):
            os.mkdir(txt_loc)
   

        eloc = f"{plot_loc}/exhumed"
        if not os.path.exists(eloc):
            os.mkdir(eloc)

        exhloc = f"{txt_loc}/exhumed"
        if not os.path.exists(exhloc):
            os.mkdir(exhloc)

        pt_files = f'{txt_loc}/PT'
        npa = len(os.listdir(pt_files))
        print("Total number of particles = ", npa)

        init = pd.read_csv(f"{txt_loc}/particles_indexes.txt", sep="\s+")
  
        exh = np.zeros(npa)-1
        exh[:] = np.nan
        pal1 = plt.get_cmap('viridis')
        norm = plt.Normalize(init["init_x"].min()/1e3, init["init_x"].max()/1e3)
        line_colors = pal1(norm(init["init_x"]/1e3))

        for i in range(0,npa):
            p = pd.read_csv(f"{pt_files}/pt_part_{i}.txt", sep="\s+")
            df = p[p["P"] <= 4.5]
            if not df.tail(2).depth.is_monotonic_decreasing:
                if (p.P > 0.1).any():
                    exh[i] = i
                #     plt.plot(p["T"] - 273, p["P"], linewidth = 0.5, c = line_colors[i])
                # plt.savefig(f"{plot_loc}/non_monotonic.png", dpi = 1000)
                    
            else:
                if (p.P > 1.5).any(): 
                    plt.plot(p["T"] - 273, p["P"], linewidth = 0.5, c = line_colors[i])
        # plt.colorbar()
        s_m = cm.ScalarMappable(cmap = pal1, norm = norm)
        s_m.set_array([])
        plt.colorbar(s_m, label = "Initial x position (km)")
        plt.ylim(0,4.5)
        plt.xlim(0,1000)
        plt.xlabel("T ($^\circ$C)")
        plt.ylabel("P (GPa)")
        plt.savefig(f"{plot_loc}/subducted.png", dpi = 1000)
        plt.close()
    
        exh = pd.Series(exh[~np.isnan(exh)])
        print("num of potentially exhumed = ", len(exh), " particles")
        ts = int(len(time_array))
        P = np.zeros(ts-1)
        T = np.zeros(ts-1)
        P[:]=np.nan
        T[:]=np.nan
        maxp = np.zeros((len(exh), 2))
        maxt = np.zeros((len(exh), 2))
        pal2 = plt.get_cmap('viridis',int(len(exh)/freq))


        threshold = .09
        filename = f"{txt_loc}/maxPT.txt"
        maxx = open(filename,"w+")
        maxx.write("part maxPP maxPT maxTP maxTT terrain lithology\n")

        count = 0
        fig, ax = plt.subplots(1, 3, figsize = (15,5))
        for ind_j, j in tqdm(enumerate(exh[::freq])):
            
            e = pd.read_csv(f"{pt_files}/pt_part_{int(j)}.txt", sep="\s+")
            e["terrain"] = 0

            for ind_c, c in enumerate(compositions):
                weight = ind_c + 1 
                e["terrain"] += weight * e[c]
            e["terrain"] = e["terrain"].astype(int)

            lithology_classes = e.terrain.unique()
            composition_mapping = {ind_c + 1: c for ind_c, c in enumerate(compositions)}
            e["lithology"] = e["terrain"].map(composition_mapping)

           
            if len(e.P) >= 10:
                if (e.depth.iat[-1] - e.depth.min()) >= 5.:
                # if e.depth.min() > (e.depth.max() - 30.):
                    # print(e.depth.max())
                    e["Plith"] = (e["y"].max()- e["y"])*9.81*3300/1e9
                    idxp = e["Plith"].idxmax()
                    idxt = e["T"].idxmax()
            
                    maxp[ind_j,0] = e["Plith"].iloc[idxp]
                    maxp[ind_j,1] = e["T"].iloc[idxp] - 273.
                    maxt[ind_j,0] = e["Plith"].iloc[idxt]
                    maxt[ind_j,1] = e["T"].iloc[idxt] - 273.

                    

                 # and abs(trench[-1] - a.x.iloc[-1]) <= 50.e3:
                    if e.P.iloc[0]<0.2:
                        count = count+1
                        maxx.write("%.0f %.3f %.3f %.3f %.3f %.0f %s\n" % (j, maxp[ind_j,0], maxp[ind_j,1], maxt[ind_j,0], maxt[ind_j,1], e["terrain"].iloc[0], e["lithology"].iloc[0]))
                        

                        ax[0].plot(e["T"] -273, e["P"], c = pal2(ind_j), linewidth = 1)
                        ax[0].scatter(e["T"].iloc[0] -273, e["P"].iloc[0], s = 10)
                        ax[0].set_xlabel("T ($^\circ$C)")
                        ax[0].set_ylabel("P (GPa)")
                        ax[0].set_ylim(0,e["P"].max() + 0.3)
                        ax[0].set_ylim(0, 2.5)
                        ax[0].set_title("Full Pressure")

                        ax[1].plot(e["T"] -273, e["Plith"], c = pal2(ind_j), linewidth = 1)
                        ax[1].set_xlabel("T ($^\circ$C)")
                        ax[1].set_ylabel("P (GPa)")
                        ax[1].set_ylim(0, 2.5)
                        ax[1].set_title("Lithostatic Pressure")
                        ax2 = ax[1].twinx()
                        ax2.set_ylim(0,76)
                        ax2.set_ylabel("depth (km)")

                        ax[2].plot((e["x"])/1.e3, (e["y"].max() - e["y"])/1.e3, c = pal2(ind_j), linewidth = 1)
                        ax[2].set_xlabel("x (km)")
                        ax[2].set_ylabel("y (km)")
                        ax[2].set_ylim(50, -5)
                        ax[2].set_title("Particle trajectory")
            
                        


                        exname = f"{exhloc}/exhumed_{int(j)}.txt"
                        expar = open(exname,"w+")
                        expar.write("id time x depth P Plith T terrain lithology\n")
                        for index, row in e.iterrows():
                            expar.write("%.0f %.3f %.3f %.3f %.3f %.3f %.3f %.0f %s\n" % (row["id"], row["time"], row["x"], row["y"], row["P"], row["Plith"], row["T"], row["terrain"], row["lithology"]))
                        expar.close()
                        
                        # count = count+1
        plt.savefig(f"{plot_loc}/potential_exhum.png", dpi = 1000)     
        plt.close()
        maxx.close()
        
        fig, ax = plt.subplots(1, 3, figsize = (15,5))
        # mycmap = colors.ListedColormap(['magenta', 'pink'])
        maxconds = pd.read_csv(f"{txt_loc}/maxPT.txt", sep="\s+")
        sns.scatterplot(ax=ax[0], data = maxconds, x = "maxPT", y = "maxPP", hue = "lithology")
        ax[0].set_xlabel("T ($^\circ$C)")
        ax[0].set_ylabel("P (GPa)")
        ax[0].set_title("Max P")

        
        sns.scatterplot(ax=ax[1], data = maxconds, x = "maxTT", y = "maxTP", hue = "lithology")
        ax[1].set_xlabel("T ($^\circ$C)")
        ax[1].set_ylabel("P (GPa)")
        ax[1].set_title("Max T")


        ax[2].pie(maxconds.lithology.value_counts(), labels = maxconds.lithology.unique(), autopct='%1.1f%%')
        ax[2].set_title("Exhumed lithology")
        fig.suptitle(f"Total number of exhumed particles = {len(maxconds)} - {(len(maxconds)/npa) * 100 :.1f} %")
        fig.tight_layout()
        plt.savefig(f"{plot_loc}/max_PT_conditions.png", dpi = 1000)
        plt.close()


        # palt = plt.get_cmap('viridis',int(len(time_array)))

        # for name in os.listdir(exhloc):
        #     nameplot = name.split(".")[0]
        #     maxconds = pd.read_csv(f"{exhloc}/{name}", sep="\s+")
        #     fig1, ax1 = plt.subplots(1, 3, figsize = (15,5), layout = 'tight')


        #     plt.suptitle(f"Lithology: {maxconds['lithology'].iloc[0]}", fontsize = 16)

        #     ax1[0].scatter(maxconds["T"]-273, maxconds["P"], c = maxconds["time"], cmap = palt)
        #     # ax1[0].annotate(maxconds[maxconds["Plith"] >=0.00]["time"].iloc[0].astype(int), (maxconds[maxconds["Plith"] >= 0.00]["T"].iloc[0]-273, maxconds[maxconds["Plith"] >= 0.1]["P"].iloc[0]))
        #     ax1[0].annotate(maxconds["time"].iloc[-1].astype(int), (maxconds["T"].iloc[-1]-273, maxconds["Plith"].iloc[-1]))
        #     for k in range(len(maxconds)):
        #         if maxconds["P"].iloc[k] == maxconds["P"].max():
        #             ax1[0].scatter(maxconds["T"].iloc[k]-273, maxconds["P"].iloc[k], c = 'black', marker="x")
        #             ax1[0].annotate(maxconds["time"].iloc[k].astype(int), (maxconds["T"].iloc[k]-273, maxconds["P"].iloc[k]))
        #         elif maxconds["T"].iloc[k] == maxconds["T"].max():
        #             ax1[0].scatter(maxconds["T"].iloc[k]-273, maxconds["P"].iloc[k], c = 'black', marker="*")
        #             ax1[0].annotate(maxconds["time"].iloc[k].astype(int), (maxconds["T"].iloc[k]-273, maxconds["P"].iloc[k]))
        #     ax1[0].set_xlabel("T ($^\circ$C)")
        #     ax1[0].set_ylabel("P (GPa)")
        #     ax1[0].set_ylim(0, maxconds["P"].max() + 0.1)

            
        
        #     sc=ax1[1].scatter(maxconds["T"]-273, maxconds["Plith"], c = maxconds["time"], cmap = palt)
        #     # ax1[1].annotate(maxconds[maxconds["Plith"] >= 0.0]["time"].iloc[0].astype(int), (maxconds[maxconds["Plith"] >= 0.0]["T"].iloc[0]-273, maxconds[maxconds["Plith"] >= 0.1]["Plith"].iloc[0]))
        #     ax1[1].annotate(maxconds["time"].iloc[-1].astype(int), (maxconds["T"].iloc[-1]-273, maxconds["Plith"].iloc[-1]))

        #     for k in range(len(maxconds)):
        #         if maxconds["Plith"].iloc[k] == maxconds["Plith"].max():
        #             ax1[1].scatter(maxconds["T"].iloc[k]-273, maxconds["Plith"].iloc[k], c = 'black', marker="x")
        #             ax1[1].annotate(maxconds["time"].iloc[k].astype(int), (maxconds["T"].iloc[k]-273, maxconds["Plith"].iloc[k]))
        #         elif maxconds["T"].iloc[k] == maxconds["T"].max():
        #             ax1[1].scatter(maxconds["T"].iloc[k]-273, maxconds["Plith"].iloc[k], c = 'black', marker="*")
        #             ax1[1].annotate(maxconds["time"].iloc[k].astype(int), (maxconds["T"].iloc[k]-273, maxconds["Plith"].iloc[k]))
        #     ax1[1].set_xlabel("T ($^\circ$C)")
        #     ax1[1].set_ylabel("$P_{lith}$ (GPa)")
        #     ax1[1].set_ylim(0, maxconds["P"].max() + 0.1)
        #     axd = ax1[1].twinx()
        #     dmax = ((maxconds["P"].max() + 0.1)*1e9)/(3300*9.81)
        #     axd.set_ylim(0,dmax/1e3)
        #     axd.set_ylabel("depth (km)")


        #     ax1[2].scatter(maxconds["x"]/1.e3, (maxconds["depth"].max() - maxconds["depth"])/1.e3, c = maxconds["time"], cmap = palt)
        #     # ax1[2].annotate(maxconds[maxconds["Plith"] >=0.1]["time"].iloc[0].astype(int), (maxconds[maxconds["Plith"] >= 0.1]["T"].iloc[0]-273, maxconds[maxconds["Plith"] >= 0.1]["P"].iloc[0]))
        #     # ax1[2].annotate(maxconds["time"].iloc[-1].astype(int), (maxconds["T"].iloc[-1]-273, maxconds["Plith"].iloc[-1]))
        #     # for k in range(len(maxconds)):
        #     #     if maxconds["P"].iloc[k] == maxconds["P"].max():
        #     #         ax1[0].scatter(maxconds["T"].iloc[k]-273, maxconds["P"].iloc[k], c = 'black', marker="x")
        #     #         ax1[0].annotate(maxconds["time"].iloc[k].astype(int), (maxconds["T"].iloc[k]-273, maxconds["P"].iloc[k]))
        #     #     elif maxconds["T"].iloc[k] == maxconds["T"].max():
        #     #         ax1[0].scatter(maxconds["T"].iloc[k]-273, maxconds["P"].iloc[k], c = 'black', marker="*")
        #     #         ax1[0].annotate(maxconds["time"].iloc[k].astype(int), (maxconds["T"].iloc[k]-273, maxconds["P"].iloc[k]))
        #     ax1[2].set_xlabel("x (km)")
        #     ax1[2].set_ylabel("y (km)")
        #     ax1[2].set_ylim((maxconds["depth"].max() - maxconds["depth"].min())/1.e3 + 2, -5)
            
        #     plt.savefig(f"{eloc}/{nameplot}.png", dpi = 1000)
        #     plt.close()
        

        

        print("num of exhumed = ", count, " particles")





if __name__ == "__main__":
    main()


