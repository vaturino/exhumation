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

def collect_initial_particles(trench: float, min_d: float, max_d: float, init):
    data = init[(init["Points:0"] < trench - min_d) & (init["Points:0"] > trench - max_d) & (init["opc"] == 0)]
    return data

def filter_particles_by_sample(model_path, sample_ids, compositions, tr):
    """
    Filters particles in the last timestep that are in the sample IDs.
    Args:
        model_path (str): Path to the model folder.
        sample_ids (set): Set of particle IDs from the sample.
        compositions (list): List of compositions to track.
        tr (float): Threshold for composition.
    Returns:
        pd.DataFrame: DataFrame of particles matching the criteria.
    """
    last_file_path = os.path.join(model_path, 'particles/full.130.gzip')
    if not os.path.exists(last_file_path):
        print(f"File not found: {last_file_path}")
        return pd.DataFrame()

    # Load the last timestep
    last_timestep_data = load_data(model_path, '/particles/full.130.gzip', compositions, tr)
    
    # Filter by IDs in sample
    filtered_data = last_timestep_data[last_timestep_data["id"].isin(sample_ids)]
    return filtered_data


def main():
    json_file = "/home/vturino/PhD/projects/exhumation/pyInput/velocity_comparison.json"
    with open(json_file) as json_file:
        file = json.load(json_file)

    csvs_loc = '/home/vturino/Vale_nas/exhumation/gz_outputs/'
    plot_loc = '/home/vturino/PhD/projects/exhumation/plots/single_models/'
    compositions = file['compositions']
    track = file['track']
    tr = 0.5   

    # Process the first model for sampling
    model = file['models'][0]
    plot_folder = file['plot_folder'][0]
    initial = load_data(csvs_loc + model, '/particles/full.1.gzip', compositions, tr)
    trench = get_trench_position_from_op(initial)
    dataset = collect_initial_particles(trench, 0.e3, 350.e3, initial)
    dataset["lithology_composition"] = [0] * len(dataset)

    for ind_c, c in enumerate(track):
        weight = ind_c + 1
        dataset["lithology_composition"] += weight * dataset[c].round()

    dataset["lithology_composition"] = dataset["lithology_composition"].round().astype(int)

    # Sample 25% of the dataset
    msk = np.random.rand(len(dataset)) <= 0.25
    sample = dataset[msk]
    sample_ids = set(sample["id"])  # Extract sample IDs for filtering
    # get the count of each class
    lithology_classes = sample.lithology_composition.unique()

    # Create lithology mapping
    lithology_mapping = {lithology_classes[k-1]: track[lithology_classes[k-1]-1] for k in range(1, len(track)+1)}

    # Map the numbers in 'lithology_composition' to their corresponding composition names
    sample['lithology'] = sample['lithology_composition'].map(lithology_mapping)


    # Loop through each model to find matching particles
    model_ids = []
    for model in file['models']:
        model_path = os.path.join(csvs_loc, model)
        matching_particles = filter_particles_by_sample(model_path, sample_ids, compositions, tr)
        model_ids.append(set(matching_particles["id"]))
        

        

 
    # Find common IDs across all models
    common_ids = set.intersection(*model_ids)

    # Generate the final filtered dataset for each model using common IDs
    results = {}
    for model in file['models']:
        model_path = os.path.join(csvs_loc, model)
        matching_particles = filter_particles_by_sample(model_path, common_ids, compositions, tr)

        # Map the 'id' to its corresponding 'lithology'
        matching_particles["lithology"] = matching_particles["id"].map(sample.set_index('id')['lithology'])
        matching_particles["lithology_composition"] = matching_particles["id"].map(sample.set_index('id')['lithology_composition'])
        #Reset index
        matching_particles.reset_index(drop=True, inplace=True)


        # Filter only the columns we care about
        fixed_columns = ['id', 'initial position:0', 'initial position:1', 'lithology_composition', 'lithology']
        compo_columns = [f"initial {c}" for c in track]
        
        matching_particles = matching_particles[fixed_columns + compo_columns]

        # Rename the columns:
        matching_particles = matching_particles.rename(columns={"initial position:0": "init_x", "initial position:1": "init_y"})
        matching_particles = matching_particles.rename(columns={f"initial {c}": f"init_{c}" for c in track})
        matching_particles = matching_particles.sort_values(by=['id']).reset_index(drop=True)
        #save index as a new column "particle"
        matching_particles.insert(0, 'particle', matching_particles.index)
        # Save the filtered data
        results[model] = matching_particles

    # Save results and visualize
    for model, data in results.items():
        print(f"Model: {model}, Common Matching Particles: {len(data)}")
        save_path = f"{plot_loc}/{model}/txt_files/particles_indexes.csv"
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        data.to_csv(save_path, index=False)

        # Scatter plot with numeric 'lithology_composition' values
        plt.scatter(data["init_x"], data["init_y"], c=data["lithology_composition"], cmap='viridis')
        plt.title(f"Particles for {model}")
        plt.xlabel("Initial X Position")
        plt.ylabel("Initial Y Position")
        plt.colorbar(label="Lithology Composition")
        plt.savefig(f"{plot_loc}/{model}/lithology_scatter.png")
        plt.clf()

if __name__ == "__main__":
    main()
