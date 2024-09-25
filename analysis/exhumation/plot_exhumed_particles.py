#! /usr/bin/python3

### Script that loops through all of the particles and identifies and plots the exhumed particles ###


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

    # Read the json file
    with open(f"{json_loc}{args.json_file}") as json_file:
            configs = json.load(json_file)
            setting = 'none'
    compositions = configs['compositions']
    # compositions = [c for c in configs['compositions'] if c != 'opc']
    
    # print(compositions)
    m = configs["models"][0]

    # Create the folders to save the plots
    plot_loc = f"/home/vturino/PhD/projects/exhumation/plots/single_models/{configs['models'][0]}"
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


     # Get the time array and the number of time steps
    time_array = np.zeros((len(os.listdir(f"{csvs_loc}{m}/particles")),2))   
    stat = pd.read_csv(f"{models_loc}{m}/statistics",skiprows=configs['head_lines'],sep='\s+',header=None)
    time_array = grab_dimTime_particles(f"{csvs_loc}{m}/particles", stat, time_array, configs['head_lines']-1)
    ts = int(len(time_array))
    init = pd.read_csv(f"{txt_loc}/particles_indexes.txt", sep="\s+")

    pal1 = plt.get_cmap('viridis')
    norm = plt.Normalize(init["init_x"].min()/1e3, init["init_x"].max()/1e3)
    line_colors = pal1(norm(init["init_x"]/1e3))

    # Create the array with the exhumed particles indexes and the colorbar for plotting
    exh = pd.DataFrame(columns=["id", "maxPP", "maxPT", "maxTP", "maxTT", "lithology"], index=range(npa))
    stag = pd.DataFrame(columns=["id", "maxPP", "maxPT", "maxTP", "maxTT", "lithology"], index=range(npa))
    pal1 = plt.get_cmap('viridis')
    norm = plt.Normalize(init["init_x"].min()/1e3, init["init_x"].max()/1e3)
    line_colors = pal1(norm(init["init_x"]/1e3))
    s_m = cm.ScalarMappable(cmap = pal1, norm = norm)
    s_m.set_array([])

    subducted = 0
    exhumed = 0
    stagnant = 0


    composition_mapping = {ind_c +1: c for ind_c, c in enumerate(compositions)}


    # Loop through all the time steps and particles to identify the exhumed particles
    f1,a1 = plt.subplots(1, 2, figsize=(15, 5))
    f2,a2 = plt.subplots(1, 1)
    f3,a3 = plt.subplots(1, 2, figsize=(15, 5))
    for p in tqdm(range(0, npa)):
        pt_single = pd.read_csv(f"{pt_files}/pt_part_{p}.txt", sep="\s+")
        pt_single["Plith"] = (pt_single["depth"].max()- pt_single["depth"])*1e3*9.81*3300/1e9
        pt_single["terrain"] = 0

        weight = 0

        for ind_c, c in enumerate(compositions):
            weight = ind_c + 1 
            pt_single["terrain"] += weight * pt_single[c].round()
        pt_single["terrain"] = pt_single["terrain"].round()
        
        

        pt_single["lithology"] = pt_single["terrain"].map(composition_mapping)
        pt_single["lithology"] = pt_single["lithology"].iloc[-1]


        if pt_single["P"].max() > 3.0:
            subducted += 1
            a2.plot(pt_single["T"]-273., pt_single["P"], color=line_colors[p])
            a2.set_xlabel("Temperature (C)")
            a2.set_ylabel("Pressure (GPa)")
            a2.set_title("Subducted particles")
            a2.set_ylim(0, 4.5)
            a2.set_xlim(0, 1100)
        else: #check that the distance between the maximum depth and the last depth is at least 5 km
            if (pt_single.depth.iat[-1] - pt_single.depth.min()) >= 5.:
                exhumed += 1

                exh["id"].iloc[p] = p
                eidxp = pt_single["Plith"].idxmax()
                eidxt = pt_single["T"].idxmax()
                exh["maxPP"].iloc[p] = pt_single["Plith"].iloc[eidxp]
                exh["maxPT"].iloc[p] = pt_single["T"].iloc[eidxp] - 273.
                exh["maxTT"].iloc[p] = pt_single["T"].iloc[eidxt] - 273.
                exh["maxTP"].iloc[p] = pt_single["Plith"].iloc[eidxt]
                exh["lithology"].iloc[p] = pt_single["lithology"].iloc[eidxp]
                
                a1[0].plot(pt_single["T"]-273., pt_single["Plith"], color=line_colors[p])
                a1[0].set_xlabel("Temperature (C)")
                a1[0].set_ylabel("Pressure (GPa)")
                a1[0].set_title("Exhumed particles")


                a1[1].plot(pt_single["x"]/1.e3, pt_single["depth"].max() - pt_single["depth"], color=line_colors[p])
                a1[1].invert_yaxis()
                a1[1].set_xlabel("Distance (km)")
                a1[1].set_ylabel("Depth (km)")
                a1[1].set_title("Exhumed particles trajectory")
                

            # stagnant particles: particles for thich the distance between the maximum depth and the last depth is less than 5 km
            elif (pt_single.depth.iat[-1] - pt_single.depth.min()) <= 5. and (pt_single.depth.iat[-1] - pt_single.depth.min()) >= 0.5:
                # print("Stagnant particle = ", p)    
                # print(pt_single["serp"].iloc[-1])
                stagnant += 1
                stag["id"].iloc[p] = p

                sidxp = pt_single["Plith"].idxmax()
                sidxt = pt_single["T"].idxmax()
                stag["maxPP"].iloc[p] = pt_single["Plith"].iloc[sidxp]
                stag["maxPT"].iloc[p] = pt_single["T"].iloc[sidxp] - 273.
                stag["maxTT"].iloc[p] = pt_single["T"].iloc[sidxt] - 273.
                stag["maxTP"].iloc[p] = pt_single["Plith"].iloc[sidxt]
                stag["lithology"].iloc[p] = pt_single["lithology"].iloc[sidxp]

                a3[0].plot(pt_single["T"]-273., pt_single["Plith"], color=line_colors[p])
                a3[0].set_xlabel("Temperature (C)")
                a3[0].set_ylabel("Pressure (GPa)")
                a3[0].set_title("Stagnant particles")

                a3[1].plot(pt_single["x"]/1.e3, pt_single["depth"].max() - pt_single["depth"], color=line_colors[p])
                a3[1].invert_yaxis()
                a3[1].set_xlabel("Distance (km)")
                a3[1].set_ylabel("Depth (km)")
                a3[1].set_title("Stagnant particles trajectory")
            else:
                stagnant += 1
                a2.plot(pt_single["T"]-273., pt_single["P"], color=line_colors[p])
                a2.set_xlabel("Temperature (C)")
                a2.set_ylabel("Pressure (GPa)")
                a2.set_title("Subducted particles")
                a2.set_ylim(0, 4.5)
                a2.set_xlim(0, 1100)

    c1 = plt.colorbar(s_m, ax=a1[1])
    c1.set_label('Initial X Position (km)')           
    c2 = plt.colorbar(s_m, ax=a2)
    c2.set_label('Initial X Position (km)')
    c3 = plt.colorbar(s_m, ax=a3[1])
    c3.set_label('Initial X Position (km)')
    f1.tight_layout()
    f3.tight_layout()
    f2.tight_layout()
    f1.savefig(f"{plot_loc}/possibly_exhumed.png")
    f3.savefig(f"{plot_loc}/stagnant.png")
    f2.savefig(f"{plot_loc}/filtered_out.png")   
    plt.close()

    exh = exh.dropna()
    stag = stag.dropna()
    print("Subducted particles = ", subducted) 
    print("Exhumed particles = ", len(exh))
    print("Stagnant particles = ", len(stag))

    #Titles 
    exh_title_str = ""
    for j in range(len(exh.lithology.value_counts())):
        lithology = exh.lithology.value_counts().index[j]
        count = exh.lithology.value_counts().values[j]
        percentage = (count / npa) * 100
        exh_title_str += f"{lithology}: count = {count}, percentage = {percentage:.1f}%\n"

    stag_title_str = ""
    for j in range(len(stag.lithology.value_counts())):
        lithology = stag.lithology.value_counts().index[j]
        count = stag.lithology.value_counts().values[j]
        percentage = (count / npa) * 100
        stag_title_str += f"{lithology}: count = {count}, percentage = {percentage:.1f}%\n"





    

    # Plot the maximum conditions of the exhumed particles
    f4, a4 = plt.subplots(1, 3, figsize = (15,6))
    sns.scatterplot(ax=a4[0], data = exh, x = "maxPT", y = "maxPP", hue = "lithology", linewidth=0.2)
    a4[0].set_xlabel("T ($^\circ$C)")
    a4[0].set_ylabel("P (GPa)")
    a4[0].set_title("Max P")
    
    sns.scatterplot(ax=a4[1], data = exh, x = "maxTT", y = "maxTP", hue = "lithology", linewidth=0.2)
    a4[1].set_xlabel("T ($^\circ$C)")
    a4[1].set_ylabel("P (GPa)")
    a4[1].set_title("Max T")

    a4[2].pie(exh.lithology.value_counts(), labels = exh.lithology.unique(), autopct='%1.1f%%')
    a4[2].set_title("Exhumed lithology")
    f4.suptitle(f"Total number of exhumed particles = {len(exh)} - {(len(exh)/npa) * 100 :.1f}%")
    f4.text(0.8, 0.1, exh_title_str, ha='center', va='center')
    f4.tight_layout()
    plt.savefig(f"{plot_loc}/max_PT_conditions.png", dpi = 1000)
    plt.close()
   

    # Plot the maximum conditions of the stagnant particles
    f5, a5 = plt.subplots(1, 3, figsize = (15,6))
    sns.scatterplot(ax=a5[0], data = stag, x = "maxPT", y = "maxPP", hue = "lithology", linewidth=0.2)
    a5[0].set_xlabel("T ($^\circ$C)")
    a5[0].set_ylabel("P (GPa)")
    a5[0].set_title("Max P")

    sns.scatterplot(ax=a5[1], data = stag, x = "maxTT", y = "maxTP", hue = "lithology", linewidth=0.2)
    a5[1].set_xlabel("T ($^\circ$C)")
    a5[1].set_ylabel("P (GPa)")
    a5[1].set_title("Max T")

    a5[2].pie(stag.lithology.value_counts(), labels = stag.lithology.unique(), autopct='%1.1f%%')
    a5[2].set_title("Stagnant lithology")
    f5.suptitle(f"Total number of stagnant particles = {len(stag)} - {(len(stag)/npa) * 100 :.1f} %")
    f5.text(0.8, 0.1, stag_title_str, ha='center', va='center')
    f5.tight_layout()
    plt.savefig(f"{plot_loc}/max_PT_conditions_stagnant.png", dpi = 1000)
    plt.close()

   

    # Save the exhumed and stagnant particles indexes
    exh.to_csv(f"{txt_loc}/exhumed_particles.txt", sep="\t", index=False)
    stag.to_csv(f"{txt_loc}/stagnant_particles.txt", sep="\t", index=False)



      


    

if __name__ == "__main__":
    main()