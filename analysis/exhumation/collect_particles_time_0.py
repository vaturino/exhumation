#! /usr/bin/python3
from matplotlib.cm import get_cmap
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'
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
from pathlib import Path
path = str(Path(Path(__file__).parent.absolute()).parent.absolute())
sys.path.insert(0, path)
from libraries.functions import *
from libraries.particles import *
from libraries.exhumation import *

def get_trench_position_from_op(p, threshold = 2.7e7):
    if {"opc"}.issubset(p.columns):
        tr =  p.loc[(p['Points:0']< threshold) & (p["opc"] > 0.3) & (p["Points:1"] >=  p["Points:1"].max() - 2.e3),'Points:0'].min()
    else:
        tr =  p.loc[(p['Points:0']< threshold) & (p['op'] > 0.3) & (p["Points:1"] >=  p["Points:1"].max() - 2.e3),'Points:0'].min()
    return tr


def main():

    parser = argparse.ArgumentParser(description= 'Script that gets some models and gives the kinematic indicators')
    parser.add_argument('json_file', help='json file with model name, time at the end of computation, folder where to save the plots, legend')
    args = parser.parse_args()

    # csvs_loc =  '/home/vturino/PhD/projects/exhumation/gz_outputs/'
    # models_loc =  '/home/vturino/PhD/projects/exhumation/raw_outputs/'

    models_loc = '/home/vturino/Vale_nas/exhumation/raw_outputs/'
    csvs_loc = '/home/vturino/Vale_nas/exhumation/gz_outputs/'

    json_loc = '/home/vturino/PhD/projects/exhumation/pyInput/'

    samples = 30
    tr = 1e-2


    with open(f"{json_loc}{args.json_file}") as json_file:
            configs = json.load(json_file)
    
    for ind_m, m in tqdm(enumerate(configs['models'])): 
        

        plot_loc = f"/home/vturino/PhD/projects/exhumation/plots/single_models/{configs['models'][ind_m]}"
        txt_loc = f'{plot_loc}/txt_files'
        if not os.path.exists(txt_loc):
            os.mkdir(txt_loc)

        part_loc = f'{csvs_loc}/{m}/particles/'
        last = len(os.listdir(f"{part_loc}")) - 1

        #set variables and compositions for the model
        tr  = 0.5
        compositions = configs['compositions']
        track = configs['track']
        

        # load data for the first and last time step
        initial = load_data(csvs_loc+m, '/particles/full.1.gzip', compositions, tr)
        final = load_data(csvs_loc+m, f'/particles/full.{last}.gzip', compositions, tr)
        trench= get_trench_position_from_op(initial)
        dataset = collect_particles(trench, 0.e3, 350.e3, initial, final)
        dataset["lithology"] = [0]*len(dataset)

        # create lithology classes
        # for c in compositions:x
        #     dataset[c][dataset[c] >= 0.5] = 1
        #     dataset[c][dataset[c] < 0.5] = 0

        for ind_c, c in enumerate(track):
            weight = ind_c + 1 
            dataset["lithology"] += weight * dataset[c].round()
            
        dataset["lithology"] = dataset["lithology"].round().astype(int)

        # sample 5% of the dataset
        msk = np.random.rand(len(dataset)) <= 0.05
        sample = dataset[msk]

        # get the count of each class
        lithology_classes = sample.lithology.unique()
        

        # # create a custom colormap
        colors_for_compositions = ['pink', 'magenta', 'purple', 'darkviolet', 'indigo', 'mediumblue', 'navy']  # Adjust the number of colors
        mycmap = colors.ListedColormap(colors_for_compositions[:len(lithology_classes)])  # Limit to the number of unique classes

        # Create lithology mapping
        lithology_mapping = {lithology_classes[k-1]: track[lithology_classes[k-1]-1] for k in range(1, len(track)+1)}

        # Map the numbers in 'lithology' to their corresponding composition names
        sample['lithology_composition'] = sample['lithology'].map(lithology_mapping)
        

        # Get sizes and labels
        sizes = sample['lithology_composition'].value_counts().values
        labels = [track[lithology_classes[i]-1] for i in range(len(lithology_classes))]
    

     
        # plot the class value counts
        fig, ax = plt.subplots(1, 2, figsize=(15, 5))
        ax[0].pie(sizes, labels=labels, autopct='%1.3f%%', colors=mycmap.colors)
        ax[0].set_title('Count of particles for each lithology')

        sns.scatterplot(data=sample, x='Points:0', y='Points:1', hue='lithology_composition', palette=mycmap.colors, ax=ax[1], size=1, linewidth=0)
        # ax[1].scatter(sample['Points:0'], sample['Points:1'], c=sample['lithology'], cmap=mycmap)
        # ax[1].legend(sample['lithology_composition'])
        ax[1].set_title('Particles distribution')
        fig.suptitle(f"Total number of particles = {len(sample)}")
        plt.tight_layout()
        plt.savefig(plot_loc + '/tracked_particles.png', dpi = 1000)
        plt.close()

        npart = len(sample)
        print("number of particles = ", npart)

        filename = f"{txt_loc}/particles_indexes.txt"
        pt = open(filename,"w+")

        if compositions == ["oc", "sed", "opc", "ecl", "serp"]:
            pt.write("particle id init_x init_y init_oc init_sed init_opc init_serp lithology\n")
            for i in tqdm(range(0,npart)):
                pt.write("%.0f %.0f %.3f %.3f %.3f %.3f %.3f %.3f %s\n" % (i, sample["id"].iloc[i], sample["initial position:0"].iloc[i], sample["initial position:1"].iloc[i], sample["oc"].iloc[i], sample["sed"].iloc[i], sample["opc"].iloc[i], sample["serp"].iloc[i], sample["lithology_composition"].iloc[i]))
        elif compositions == ["oc", "sed", "opc", "ecl"]:
            pt.write("particle id init_x init_y init_oc init_sed init_opc lithology\n")
            for i in tqdm(range(0,npart)):
                pt.write("%.0f %.0f %.3f %.3f %.3f %.3f %.3f %s\n" % (i, sample["id"].iloc[i], sample["initial position:0"].iloc[i], sample["initial position:1"].iloc[i], sample["oc"].iloc[i], sample["sed"].iloc[i], sample["opc"].iloc[i], sample["lithology_composition"].iloc[i]))
        elif compositions == ['oc', 'sed', 'opc', 'ecl', 'gabbro']:
            pt.write("particle id init_x init_y init_oc init_sed init_gabbro init_opc lithology\n")
            for i in tqdm(range(0,npart)):
                pt.write("%.0f %.0f %.3f %.3f %.3f %.3f %.3f %.3f %s\n" % (i, sample["id"].iloc[i], sample["initial position:0"].iloc[i], sample["initial position:1"].iloc[i], sample["oc"].iloc[i], sample["sed"].iloc[i], sample["gabbro"].iloc[i], sample["opc"].iloc[i], sample["lithology_composition"].iloc[i]))
        elif compositions == ['oc', 'sed', 'opc', 'ecl', 'gabbro', 'serp']:
            pt.write("particle id init_x init_y init_oc init_sed init_gabbro init_opc init_serp lithology\n")
            for i in tqdm(range(0,npart)):
                pt.write("%.0f %.0f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %s\n" % (i, sample["id"].iloc[i], sample["initial position:0"].iloc[i], sample["initial position:1"].iloc[i], sample["oc"].iloc[i], sample["sed"].iloc[i], sample["gabbro"].iloc[i], sample["opc"].iloc[i], sample["serp"].iloc[i], sample["lithology_composition"].iloc[i]))
        pt.close()

if __name__ == "__main__":
    main()


