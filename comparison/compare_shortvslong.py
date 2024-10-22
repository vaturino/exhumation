#! /usr/bin/python3

# Script that compares the models contained in the json file "models.json"
# and plots the results of the comparison: the percentage of exhumed and stagnant particles
# for each model and their weighted average peak P values with their deviation.

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import json
import os
import seaborn as sns
import textwrap

def main():
    # Load json file
    with open('shortvslong.json') as f:
        models = json.load(f)

    plot_loc = "/home/vturino/PhD/projects/exhumation/plots/single_models"

    model_names = models["models_names"]
            
    exhumation_percentages = []
    stagnation_percentages = []
    lithologies_exhumed = []  
    lithologies_stagnant = []  
    exhumation_by_lithology_all_models = []
    stagnation_by_lithology_all_models = []

    # combined_data_exhumed = pd.DataFrame(columns=["model", "total"])
    # combined_data_stagnant = pd.DataFrame(columns=["model", "total"])


    for ind_m, m in enumerate(models["models"]):
        text_loc = f"{plot_loc}/{m}/txt_files"
        initial = pd.read_csv(f"{text_loc}/particles_indexes.txt", sep="\t")
        exhumed = pd.read_csv(f"{text_loc}/exhumed_particles.txt", sep="\t")
        stagnant = pd.read_csv(f"{text_loc}/stagnant_particles.txt", sep="\t")

        # Calculate percentage of exhumed and stagnant particles
        exhumed_percentage = len(exhumed) / len(initial) * 100
        stagnant_percentage = len(stagnant) / len(initial) * 100

        exhumation_percentages.append(exhumed_percentage)
        stagnation_percentages.append(stagnant_percentage)

        # Group by lithology and calculate percentages for exhumed particles
        exhumed_by_lithology = exhumed.groupby("lithology").size() / len(initial) * 100
        exhumation_by_lithology_all_models.append(exhumed_by_lithology)
        lithologies_exhumed.append(set(exhumed["lithology"].unique()))  # Store lithologies for exhumed

        # Group by lithology and calculate percentages for stagnant particles
        stagnant_by_lithology = stagnant.groupby("lithology").size() / len(initial) * 100
        stagnation_by_lithology_all_models.append(stagnant_by_lithology)
        lithologies_stagnant.append(set(stagnant["lithology"].unique()))  # Store lithologies for stagnant

    # Combine all data into a single DataFrame for exhumed and stagnant particles
    combined_data_exhumed = pd.concat(exhumation_by_lithology_all_models, axis=1).T
    combined_data_stagnant = pd.concat(stagnation_by_lithology_all_models, axis=1).T
    

    combined_data_exhumed["total"] = exhumation_percentages
    combined_data_stagnant["total"] = stagnation_percentages
    combined_data_exhumed["model"] = model_names[:len(combined_data_exhumed)]
    combined_data_stagnant["model"] = model_names[:len(combined_data_stagnant)]

    combined_data_exhumed["perctoref"] = combined_data_exhumed["total"] / combined_data_exhumed["total"].iloc[0] * 100
    combined_data_stagnant["perctoref"] = combined_data_stagnant["total"] / combined_data_stagnant["total"].iloc[0] * 100

    combined_data_exhumed["subducted"] = 100. - combined_data_exhumed["total"]
    combined_data_stagnant["subducted"] = 100. - combined_data_stagnant["total"]

    combined_data_exhumed["lost"] = combined_data_exhumed["perctoref"].iloc[0] - combined_data_exhumed["perctoref"]
    combined_data_stagnant["lost"] = combined_data_stagnant["perctoref"].iloc[0] - combined_data_stagnant["perctoref"]

    print(combined_data_exhumed)
    print(combined_data_stagnant)
  
    # Define the layout with mosaic
    f1, a1 = plt.subplot_mosaic([['upper left', 'right'],
                                ['lower left', 'lower center']],
                                figsize=(10, 5), layout="constrained")

    # Pie chart for Exhumed vs Subducted particles
    labels_exh = ['Exhumed', 'Subducted']
    sizes_exh = [combined_data_exhumed['total'].iloc[0], combined_data_exhumed['subducted'].iloc[0]]
    colors_exh = ['cornflowerblue', 'powderblue']
    a1['upper left'].pie(sizes_exh, labels=labels_exh, colors=colors_exh, autopct='%1.1f%%', startangle=90, radius = 1.1) 
    a1['upper left'].text(1.2, 0, f"Exhumed lithologies:\n", fontsize=10)
    for i, lith in enumerate(lithologies_exhumed[0]):
        a1['upper left'].text(1.2, 0 - i*0.25, f"{lith}: {combined_data_exhumed[lith].iloc[0]:.2f}%", fontsize=10)
    a1['upper left'].set_title('Dynamic slowdown:')

    # Pie chart for Stagnant vs Subducted particles
    labels_stag = ['Stagnant', 'Subducted']
    sizes_stag = [combined_data_stagnant['total'].iloc[0], combined_data_stagnant['subducted'].iloc[0]]
    colors_stag = ['olivedrab', 'powderblue']
    a1['lower left'].pie(sizes_stag, labels=labels_stag, colors=colors_stag, autopct='%1.1f%%', startangle=90, radius = 1.1)
    a1['lower left'].text(1.2, 0, f"Stagnant lithologies:\n", fontsize=10)
    for i, lith in enumerate(lithologies_exhumed[0]):
        a1['lower left'].text(1.2, 0 - i*0.25, f"{lith}: {combined_data_stagnant[lith].iloc[0]:.2f}%", fontsize=10)

    # Line plot for particle loss
    a1['right'].plot(combined_data_exhumed["model"].iloc[1:], combined_data_exhumed["lost"].iloc[1:], 
                    label="Exhumed", marker="o", c="cornflowerblue")
    a1['right'].plot(combined_data_stagnant["model"].iloc[1:], combined_data_stagnant["lost"].iloc[1:], 
                    label="Stagnant", marker="o", c="olivedrab")
    a1['right'].legend()
    a1['right'].set_ylabel("Particle loss (%)")
    a1['right'].set_xlabel("Model")
    a1['right'].set_title("Percentage of particle loss for each model")

    # Overall title
    f1.suptitle("Comparison to dynamic slowdown case", fontsize=16)
    f1.tight_layout()

    # Save the figure
    plt.savefig("velocity/particle_loss.png", dpi=300)

        


if __name__ == "__main__":
    main()
