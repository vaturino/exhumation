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
    exhumation_numbers = []
    stagnation_numbers = []
    lithologies_exhumed = []  
    lithologies_stagnant = []  
    exhumation_by_lithology_all_models = []
    stagnation_by_lithology_all_models = []
    nparts = [] 


    
    model_names = models["models_names"]

    for ind_m, m in enumerate(models["models"]):
        text_loc = f"{plot_loc}/{m}/txt_files"
        initial = pd.read_csv(f"{text_loc}/particles_indexes.txt", sep="\t")
        exhumed = pd.read_csv(f"{text_loc}/exhumed_particles.txt", sep="\t")
        stagnant = pd.read_csv(f"{text_loc}/stagnant_particles.txt", sep="\t")

        nparts.append(len(initial))

        # Calculate percentage of exhumed and stagnant particles
        exhumed_percentage = len(exhumed) / len(initial) * 100
        stagnant_percentage = len(stagnant) / len(initial) * 100

        exhumation_percentages.append(exhumed_percentage)
        stagnation_percentages.append(stagnant_percentage)

        exhumation_number = len(exhumed)
        stagnation_number = len(stagnant)


        exhumation_numbers.append(exhumation_number)
        stagnation_numbers.append(stagnation_number)

        # Group by lithology and calculate number of particles for exhumed lithologies
        exhumed_by_lithology = exhumed.groupby("lithology").size() 
        exhumation_by_lithology_all_models.append(exhumed_by_lithology)
        lithologies_exhumed.append(set(exhumed["lithology"].unique()))  # Store lithologies for exhumed


        # Group by lithology and calculate number of particles for stagnant lithologies
        stagnant_by_lithology = stagnant.groupby("lithology").size()
        stagnation_by_lithology_all_models.append(stagnant_by_lithology)
        lithologies_stagnant.append(set(stagnant["lithology"].unique()))  # Store lithologies for stagnant

        
    
    # Combine all data into a single DataFrame for exhumed and stagnant particles
    combined_data_exhumed = pd.concat(exhumation_by_lithology_all_models, axis=1).T
    combined_data_stagnant = pd.concat(stagnation_by_lithology_all_models, axis=1).T

    combined_data_exhumed["initial"] = nparts
    combined_data_stagnant["initial"] = nparts
    
    combined_data_exhumed["values"] = models[f"velocity_values"]
    combined_data_stagnant["values"] = models[f"velocity_values"]

    combined_data_exhumed["total_perc"] = exhumation_percentages
    combined_data_stagnant["total_perc"] = stagnation_percentages

    combined_data_exhumed["total_abs"] = exhumation_numbers
    combined_data_stagnant["total_abs"] =  stagnation_numbers

    combined_data_exhumed["model"] = model_names[:len(combined_data_exhumed)]
    combined_data_stagnant["model"] = model_names[:len(combined_data_stagnant)]


    combined_data_exhumed["perctoref"] = (combined_data_exhumed["total_abs"] - combined_data_exhumed["total_abs"].iloc[0]) / combined_data_exhumed["total_abs"].iloc[0] * 100
    combined_data_stagnant["perctoref"] = (combined_data_stagnant["total_abs"] - combined_data_stagnant["total_abs"].iloc[0]) / combined_data_stagnant["total_abs"].iloc[0] * 100

    combined_data_exhumed["subducted"] = 100 - combined_data_exhumed["total_perc"] 
    combined_data_stagnant["subducted"] = 100 - combined_data_stagnant["total_perc"]
    combined_data_exhumed["overall_subducted"] = 100 - combined_data_exhumed["total_perc"] - combined_data_stagnant["total_perc"]

        

    suffix = "_perc"
    for lith in lithologies_exhumed[0]:
            combined_data_exhumed[lith + suffix] = (combined_data_exhumed[lith] - combined_data_exhumed[lith].iloc[0]) / combined_data_exhumed[lith].iloc[0] * 100
    for lith in lithologies_stagnant[0]:
        combined_data_stagnant[lith + suffix] = (combined_data_stagnant[lith] - combined_data_stagnant[lith].iloc[0]) / combined_data_stagnant[lith].iloc[0] * 100


   

    # Pie chart for Exhumed vs Subducted particles
    f1, a1 = plt.subplots(1, 2, figsize=(12, 5), gridspec_kw={'width_ratios': [1, 2]})
    
    
    # Pie chart for Exhumed vs Subducted particles
    labels = ['Exhumed', 'Subducted', 'Stagnant']
    sizes_overall = [combined_data_exhumed['total_perc'].iloc[0], combined_data_exhumed['overall_subducted'].iloc[0], combined_data_stagnant['total_perc'].iloc[0]]
    colors_overall = ['cornflowerblue', 'powderblue', 'olivedrab']

    #wrap names for better visualization
    model_names = [textwrap.fill(name, width=15) for name in model_names]
    
    init = combined_data_exhumed["initial"].iloc[0]


    a1[0].pie(sizes_overall, labels=labels, colors=colors_overall,
                autopct='%1.1f%%', startangle=90, radius=1, pctdistance=0.85)
    
    a1[0].text(-1.7, 0.15, "Exhumed particles:\n", fontsize=10, ha='center', va='center', weight='bold')
    a1[0].set_title(f"Percentages for {model_names[0]} model", fontsize=12)
    for i, lith in enumerate(lithologies_exhumed[0]):
        a1[0].text(-1.7, 0 - i*0.15, f"{lith}: {combined_data_exhumed[lith].iloc[0]/init * 100:.2f}%", fontsize=10, ha='center', va='center')
    a1[0].text(1.7, 0.15, "Stagnant particles:\n", fontsize=10, ha='center', va='center', weight='bold')
    for i, lith in enumerate(lithologies_stagnant[0]):
        a1[0].text(1.7, 0 - i*0.15, f"{lith}: {combined_data_stagnant[lith].iloc[0]/init * 100:.2f}%", fontsize=10, ha='center', va='center')


    # Line plot for particle loss: on x axis the value of "velocity_values" in the json file and on y axis the percentage of particle loss
    a1[1].plot(combined_data_exhumed["values"].iloc[1:], combined_data_exhumed["perctoref"].iloc[1:], 
                    label="Exhumed", marker="o", c="cornflowerblue")
    a1[1].plot(combined_data_stagnant["values"].iloc[1:], combined_data_stagnant["perctoref"].iloc[1:],
                    label="Stagnant", marker="o", c="olivedrab")
    a1[1].axhline(0, color='grey', lw=1, linestyle = "--")
    a1[1].legend()
    a1[1].set_xticks(combined_data_exhumed["values"].iloc[1:])
    a1[1].set_xticklabels(combined_data_exhumed["model"].iloc[1:])
    a1[1].set_ylabel("Particle loss (%)")
    a1[1].set_xlabel("Percentage of steady state velocity")
    a1[1].set_title("Percentage of particle loss for each model")
    a1[1].invert_xaxis()
    




    # Overall title
    f1.suptitle("Comparison to dynamic slowdown case", fontsize=16)
    f1.tight_layout()

    # Save the figure
    plt.savefig("plots/velocity/particle_loss.png", dpi=300)

        


if __name__ == "__main__":
    main()
