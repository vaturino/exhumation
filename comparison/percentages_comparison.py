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

    colors_tin = {
        "sed": "midnightblue",
        "oc": "saddlebrown",
        "ecl": "darkgreen",
        "serp": "maroon"
    }

    # Load json file
    with open('models.json') as f:
        models = json.load(f)

    plot_loc = "/home/vturino/PhD/projects/exhumation/plots/single_models"

    tests = ["velocity"] #, "viscosity", "friction", "serpentinization"]
    names = ["velocity_names"] #, "viscosity_names", "friction_names", "serpentinization_names"]
            
    exhumation_percentages = []
    stagnation_percentages = []
    lithologies_exhumed = []  
    lithologies_stagnant = []  
    exhumation_by_lithology_all_models = []
    stagnation_by_lithology_all_models = []

    combined_data_exhumed = pd.DataFrame(columns=["model", "total"])
    combined_data_stagnant = pd.DataFrame(columns=["model", "total"])


    for idx, test in enumerate(tests):

        if test in models:
            model_names = models[names[tests.index(test)]]
            all_data_exhumed = []  # List to store all dataframes for exhumed particles
            all_data_stagnant = []  # List to store all dataframes for stagnant particles

            for ind_m, m in enumerate(models[test]):
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

        combined_data_exhumed["pervel"] = models["vel_percentages"]
        combined_data_stagnant["pervel"] = models["vel_percentages"]

        combined_data_exhumed["total"] = exhumation_percentages
        combined_data_stagnant["total"] = stagnation_percentages
        combined_data_exhumed["model"] = model_names[:len(combined_data_exhumed)]
        combined_data_stagnant["model"] = model_names[:len(combined_data_stagnant)]

        combined_data_exhumed["perctoref"] = (combined_data_exhumed["total"] - combined_data_exhumed["total"].iloc[0]) / combined_data_exhumed["total"].iloc[0] * 100
        combined_data_stagnant["perctoref"] = (combined_data_stagnant["total"] - combined_data_stagnant["total"].iloc[0]) / combined_data_stagnant["total"].iloc[0] * 100

        combined_data_exhumed["subducted"] = 100 - combined_data_exhumed["total"]
        combined_data_stagnant["subducted"] = 100 - combined_data_stagnant["total"]

        

        

        suffix = "_perc"
        for lith in lithologies_exhumed[0]:
            combined_data_exhumed[lith + suffix] = (combined_data_exhumed[lith] - combined_data_exhumed[lith].iloc[0]) / combined_data_exhumed[lith].iloc[0] * 100
        for lith in lithologies_stagnant[0]:
            combined_data_stagnant[lith + suffix] = (combined_data_stagnant[lith] - combined_data_stagnant[lith].iloc[0]) / combined_data_stagnant[lith].iloc[0] * 100
        
            




        # Plot the results
        f1, a1 = plt.subplots(2, 2, figsize=(15, 10))
        
        # Pie chart for Exhumed vs Subducted particles
        labels = ['Exhumed', 'Subducted']
        sizes_exh = [combined_data_exhumed['total'].iloc[0], combined_data_exhumed['subducted'].iloc[0]]
        colors_exh = ['cornflowerblue', 'powderblue']
        sizes_stag = [combined_data_stagnant['total'][0], combined_data_stagnant['subducted'][0]]
        colors_stag = ['olivedrab', 'powderblue']


        ax_left = a1[0, 0].inset_axes([-0.1, 0.5, 0.35, 0.35])  # Left inset
        ax_right = a1[0, 0].inset_axes([0.75, 0.5, 0.35, 0.35])  # Right inset

        ax_left.pie(sizes_exh, labels=labels, colors=colors_exh,
                    autopct='%1.1f%%', startangle=90, radius=2)
        ax_left.text(0, -3, "Exhumed lithologies:\n", fontsize=10, ha='center', va='center')
        for i, lith in enumerate(lithologies_exhumed[0]):
            ax_left.text(0, -3.2 - i*0.25, f"{lith}: {combined_data_exhumed[lith].iloc[0]:.2f}%", fontsize=10, ha='center', va='center')

        
        ax_right.pie(sizes_stag, labels=labels, colors=colors_stag,
                    autopct='%1.1f%%', startangle=90, radius=2)
        ax_right.text(0, -3, "Stagnant lithologies:\n", fontsize=10, ha='center', va='center')
        for i, lith in enumerate(lithologies_stagnant[0]):
            ax_right.text(0, -3.2 - i*0.25, f"{lith}: {combined_data_stagnant[lith].iloc[0]:.2f}%", fontsize=10, ha='center', va='center')
        
        a1[0, 0].set_frame_on(False)
        a1[0, 0].tick_params(left=False, bottom=False)
        # Remove numbers on axes
        a1[0, 0].set_xticks([])
        a1[0, 0].set_yticks([])
        a1[0, 0].set_title("Constant velocity", fontsize=14)
        a1[0, 0].set_aspect('equal')


        # Plot the percentage of exhumed and stagnantparticles with respect to reference
        a1[0,1].plot(combined_data_exhumed["pervel"], combined_data_exhumed["perctoref"], marker='o', color='cornflowerblue')
        a1[0,1].plot(combined_data_stagnant["pervel"], combined_data_stagnant["perctoref"], marker='o', color='olivedrab')
        a1[0,1].axhline(y=0, color='grey', linestyle='--', zorder = 1)
        a1[0,1].set_ylabel("Percentage of particles (%)", fontsize=12)
        a1[0,1].set_xlabel("Steady state velocity percentage (%)", fontsize=12)
        a1[0,1].invert_xaxis()
        a1[0,1].set_xticks(combined_data_exhumed["pervel"])
        a1[0,1].set_title(f"Percentage of recovered particles with respect to {model_names[0]} model", fontsize=13)
        a1[0,1].legend(["Exhumed", "Stagnant"], fontsize=12)

        #plot the exhumed particles with respect to reference by lithology
        for i, lith in enumerate(lithologies_exhumed[0]):
            a1[1,0].plot(combined_data_exhumed["pervel"], combined_data_exhumed[lith+suffix], marker='o', c=colors_tin[lith])
        a1[1,0].axhline(y=0, color='grey', linestyle='--', zorder = 1)
        a1[1,0].set_ylabel("Percentage of particles (%)", fontsize=12)
        a1[1,0].set_xlabel("Steady state velocity percentage (%)", fontsize=12)
        a1[1,0].invert_xaxis()
        a1[1,0].set_xticks(combined_data_exhumed["pervel"])
        a1[1,0].set_title(f"Percentage of exhumed particles with respect to {model_names[0]} model", fontsize=13)
        a1[1,0].legend(lithologies_exhumed[0], fontsize=12)

        #plot the stagnant particles with respect to reference by lithology
        for i, lith in enumerate(lithologies_stagnant[0]):
            a1[1,1].plot(combined_data_stagnant["pervel"], combined_data_stagnant[lith+suffix], marker='o', c=colors_tin[lith])
        a1[1,1].axhline(y=0, color='grey', linestyle='--', zorder = 1)
        a1[1,1].set_ylabel("Percentage of particles (%)", fontsize=12)
        a1[1,1].set_xlabel("Steady state velocity percentage (%)", fontsize=12)
        a1[1,1].invert_xaxis()
        a1[1,1].set_xticks(combined_data_stagnant["pervel"])
        a1[1,1].set_title(f"Percentage of recovered particles with respect to {model_names[0]} model", fontsize=13)
        a1[1,1].legend(lithologies_stagnant[0], fontsize=12)


        # Overall title
        f1.suptitle(f"Comparison to {model_names[0]}", fontsize=16)
        f1.tight_layout()

        # Save the figure
        plt.savefig(f"{test}/{test}_comparison_to_ref.png", dpi=300)


if __name__ == "__main__":
    main()
