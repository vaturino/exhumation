#! /usr/bin/python3

# script that compares the models contained in the json file "models.json"
# and plots the results of the comparison: the percentage of exhumed and stagnant particles
# for each model and their avg peak P values with their deviation

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import json
import os
import seaborn as sns
import textwrap

def main():
    # Load json file
    with open('models.json') as f:
        models = json.load(f)

    plot_loc = "/home/vturino/PhD/projects/exhumation/plots/single_models"

    tests = ["velocity", "viscosity", "friction", "serpentinization"]
    names = ["velocity_names", "viscosity_names", "friction_names", "serpentinization_names"]

    for test in tests:
        if test in models:
            model_names = models[names[tests.index(test)]]
            
            exhumation_percentages = []
            stagnation_percentages = []
            lithologies_exhumed = []  # List to hold lithologies for exhumed particles
            lithologies_stagnant = []  # List to hold lithologies for stagnant particles

            exhumation_by_lithology_all_models = []
            stagnation_by_lithology_all_models = []

            avg_maxPP_exhumed = []
            std_maxPP_exhumed = []
            avg_maxPP_stagnant = []
            std_maxPP_stagnant = []

            # Store avg and std for each lithology
            lithology_maxPP_exhumed = {}
            lithology_maxPP_stagnant = {}

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

                # Calculate average peak P conditions for exhumed and stagnant particles
                avg_maxPP_exhumed.append(exhumed['maxPP'].mean())
                std_maxPP_exhumed.append(exhumed['maxPP'].std())
                avg_maxPP_stagnant.append(stagnant['maxPP'].mean())
                std_maxPP_stagnant.append(stagnant['maxPP'].std())

                # Store maxPP values for each lithology for exhumed and stagnant
                for lithology in exhumed["lithology"].unique():
                    if lithology not in lithology_maxPP_exhumed:
                        lithology_maxPP_exhumed[lithology] = []
                    lithology_maxPP_exhumed[lithology].extend(exhumed[exhumed["lithology"] == lithology]["maxPP"])

                for lithology in stagnant["lithology"].unique():
                    if lithology not in lithology_maxPP_stagnant:
                        lithology_maxPP_stagnant[lithology] = []
                    lithology_maxPP_stagnant[lithology].extend(stagnant[stagnant["lithology"] == lithology]["maxPP"])

            # Collect unique lithologies across all models
            lithologies_exhumed_all = sorted(set().union(*lithologies_exhumed))
            lithologies_stagnant_all = sorted(set().union(*lithologies_stagnant))

            # Create a color palette for lithologies
            colors_exhumed = sns.color_palette("colorblind", len(lithologies_exhumed_all))
            colors_stagnant = sns.color_palette("colorblind", len(lithologies_stagnant_all))
            lithology_colors_exhumed = dict(zip(lithologies_exhumed_all, colors_exhumed))  # Map lithologies to their colors
            lithology_colors_stagnant = dict(zip(lithologies_stagnant_all, colors_stagnant))  # Map lithologies to their colors

            # Set up the plot
            fig, ax = plt.subplots(2, 3, figsize=(15, 10))
            fig.suptitle(f"Test: {test}")

            # Bar plot for percentage of exhumed and stagnant particles by model
            x = np.arange(len(model_names))  # X-axis locations for the models
            width = 0.35  # Width of the bars

            ax[0, 0].bar(x - width/2, exhumation_percentages, width, label='Exhumation')
            ax[0, 0].bar(x + width/2, stagnation_percentages, width, label='Stagnation')

            # Configure ax[0, 0] layout
            wrapped_model_names = [textwrap.fill(name, 10) for name in model_names]  # Wrap model names
            ax[0, 0].set_xticks(x)
            ax[0, 0].set_xticklabels(wrapped_model_names, ha="center")
            ax[0, 0].set_ylabel("%")
            ax[0, 0].set_title("Exhumation and Stagnation by Model")
            ax[0, 0].legend()

            # Bar plot for exhumation by lithology for each model (ax[0, 1])
            lithology_width = width / len(lithologies_exhumed_all)  # Adjusted width for lithology bars
            for i, exhumed_by_lithology in enumerate(exhumation_by_lithology_all_models):
                for j, lithology in enumerate(lithologies_exhumed_all):
                    lithology_percentage = exhumed_by_lithology.get(lithology, 0)
                    ax[0, 1].bar(x[i] + j * lithology_width, lithology_percentage, lithology_width, 
                                  color=lithology_colors_exhumed[lithology], label=lithology if i == 0 else "")

            ax[0, 1].set_xticks(x)
            ax[0, 1].set_xticklabels(wrapped_model_names, ha="center")
            ax[0, 1].set_ylabel("%")
            ax[0, 1].set_title("Exhumation by Lithology")
            ax[0, 1].legend(title="Lithology")

            # Bar plot for stagnation by lithology for each model (ax[0, 2])
            lithology_width_stagnant = width / len(lithologies_stagnant_all)  # Adjusted width for stagnant lithology bars
            for i, stagnant_by_lithology in enumerate(stagnation_by_lithology_all_models):
                for j, lithology in enumerate(lithologies_stagnant_all):
                    lithology_percentage = stagnant_by_lithology.get(lithology, 0)
                    ax[0, 2].bar(x[i] + j * lithology_width_stagnant, lithology_percentage, lithology_width_stagnant, 
                                  color=lithology_colors_stagnant[lithology], label=lithology if i == 0 else "")

            ax[0, 2].set_xticks(x)
            ax[0, 2].set_xticklabels(wrapped_model_names, ha="center")
            ax[0, 2].set_ylabel("%")
            ax[0, 2].set_title("Stagnation by Lithology")
            ax[0, 2].legend(title="Lithology")

            # Scatter plot for average peak P conditions (ax[1, 0])
            ax[1, 0].errorbar(x, avg_maxPP_exhumed, yerr=std_maxPP_exhumed, fmt='o', label='Exhumed', color='blue', capsize=5)
            ax[1, 0].errorbar(x, avg_maxPP_stagnant, yerr=std_maxPP_stagnant, fmt='o', label='Stagnant', color='orange', capsize=5)

            ax[1, 0].set_xticks(x)
            ax[1, 0].set_xticklabels(wrapped_model_names, ha="center")
            ax[1, 0].set_ylabel("Average Peak P")
            ax[1, 0].set_title("Average Peak P Conditions")
            ax[1, 0].legend()






            ############################################################################################
            ### Scatter plot for average maxPP grouped by lithology (ax[1, 1]) for exhumed particles ###
            ############################################################################################
            avg_maxPP_by_lithology = {}
            std_maxPP_by_lithology = {}

            # Initialize dictionaries to hold maxPP values grouped by lithology
            for lithology in lithologies_exhumed_all:
                avg_maxPP_by_lithology[lithology] = []
                std_maxPP_by_lithology[lithology] = []

            # Calculate average and std for each lithology
            for ind_m, m in enumerate(models[test]):
                text_loc = f"{plot_loc}/{m}/txt_files"
                exhumed = pd.read_csv(f"{text_loc}/exhumed_particles.txt", sep="\t")

                # Group by lithology
                grouped = exhumed.groupby("lithology")['maxPP'].agg(['mean', 'std']).reset_index()
                
                # Store the average and std values for each lithology
                for index, row in grouped.iterrows():
                    lithology = row['lithology']
                    avg_maxPP_by_lithology[lithology].append(row['mean'])
                    std_maxPP_by_lithology[lithology].append(row['std'])
                
                # Append NaN for missing lithologies in other models
                for lithology in lithologies_exhumed_all:
                    if lithology not in grouped['lithology'].values:
                        avg_maxPP_by_lithology[lithology].append(np.nan)
                        std_maxPP_by_lithology[lithology].append(np.nan)

            # Prepare data for plotting
            for j, lithology in enumerate(lithologies_exhumed_all):
                ax[1, 1].errorbar(x, avg_maxPP_by_lithology[lithology], 
                                  yerr=std_maxPP_by_lithology[lithology], 
                                  fmt='o', label=lithology, capsize=5)

            ax[1, 1].set_xticks(x)
            ax[1, 1].set_xticklabels(wrapped_model_names, ha="center")
            ax[1, 1].set_ylabel("Average MaxPP")
            ax[1, 1].set_title("Average MaxPP by Lithology")
            ax[1, 1].legend(title="Lithology")




            ############################################################################################
            ### Scatter plot for average maxPP grouped by lithology (ax[1, 2]) for stagnant particles ###
            ############################################################################################
            avg_maxPP_by_lithology = {}
            std_maxPP_by_lithology = {}

            # Initialize dictionaries to hold maxPP values grouped by lithology
            for lithology in lithologies_stagnant_all:
                avg_maxPP_by_lithology[lithology] = []
                std_maxPP_by_lithology[lithology] = []

            # Calculate average and std for each lithology
            for ind_m, m in enumerate(models[test]):
                text_loc = f"{plot_loc}/{m}/txt_files"
                stagnant = pd.read_csv(f"{text_loc}/stagnant_particles.txt", sep="\t")

                # Group by lithology
                grouped = stagnant.groupby("lithology")['maxPP'].agg(['mean', 'std']).reset_index()
                
                # Store the average and std values for each lithology
                for index, row in grouped.iterrows():
                    lithology = row['lithology']
                    avg_maxPP_by_lithology[lithology].append(row['mean'])
                    std_maxPP_by_lithology[lithology].append(row['std'])
                
                # Append NaN for missing lithologies in other models
                for lithology in lithologies_stagnant_all:
                    if lithology not in grouped['lithology'].values:
                        avg_maxPP_by_lithology[lithology].append(np.nan)
                        std_maxPP_by_lithology[lithology].append(np.nan)

            # Prepare data for plotting
            for j, lithology in enumerate(lithologies_stagnant_all):
                ax[1, 2].errorbar(x, avg_maxPP_by_lithology[lithology], 
                                  yerr=std_maxPP_by_lithology[lithology], 
                                  fmt='o', label=lithology, capsize=5)
                
            ax[1, 2].set_xticks(x)
            ax[1, 2].set_xticklabels(wrapped_model_names, ha="center")
            ax[1, 2].set_ylabel("Average MaxPP")
            ax[1, 2].set_title("Average MaxPP by Lithology")
            ax[1, 2].legend(title="Lithology")
            


            # Save the plot
            plt.tight_layout()
            plt.savefig(f"combined/avg_std/{test}_comparison.png")
            plt.close()

if __name__ == "__main__":
    main()

