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

            weighted_avg_maxPP_exhumed = []
            weighted_avg_maxPP_stagnant = []

            # Store maxPP values for each lithology for exhumed and stagnant
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

                # Calculate weighted average peak P conditions for exhumed and stagnant particles
                weights_exhumed = np.ones(len(exhumed))  # Placeholder for actual weights if available
                weighted_avg_exhumed = np.average(exhumed['maxPP'], weights=weights_exhumed)
                weighted_avg_maxPP_exhumed.append(weighted_avg_exhumed)

                weights_stagnant = np.ones(len(stagnant))  # Placeholder for actual weights if available
                weighted_avg_stagnant = np.average(stagnant['maxPP'], weights=weights_stagnant)
                weighted_avg_maxPP_stagnant.append(weighted_avg_stagnant)

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

           # Prepare for scatter plot for weighted average peak P conditions (ax[1, 0])
            yerr_exhumed = []
            yerr_stagnant = []

            # Ensure x is same length as the averages for exhumed and stagnant
            x_weighted = np.arange(len(weighted_avg_maxPP_exhumed))

            # Calculate error bars (min and max) for each weighted average
            for ind_m, m in enumerate(models[test]):
                text_loc = f"{plot_loc}/{m}/txt_files"
                exhumed = pd.read_csv(f"{text_loc}/exhumed_particles.txt", sep="\t")
                stagnant = pd.read_csv(f"{text_loc}/stagnant_particles.txt", sep="\t")

                min_maxPP_exhumed = np.min(exhumed['maxPP'])
                max_maxPP_exhumed = np.max(exhumed['maxPP'])
                yerr_exhumed.append([weighted_avg_maxPP_exhumed[ind_m] - min_maxPP_exhumed, max_maxPP_exhumed - weighted_avg_maxPP_exhumed[ind_m]])

                min_maxPP_stagnant = np.min(stagnant['maxPP'])
                max_maxPP_stagnant = np.max(stagnant['maxPP'])
                yerr_stagnant.append([weighted_avg_maxPP_stagnant[ind_m] - min_maxPP_stagnant, max_maxPP_stagnant - weighted_avg_maxPP_stagnant[ind_m]])

            # Define offsets for the points
            offset = 0.1  # Adjust this value as needed

            # Plot weighted averages with error bars
            ax[1, 0].errorbar(x_weighted - offset, weighted_avg_maxPP_exhumed, 
                            yerr=np.array(yerr_exhumed).T, 
                            fmt='o', label='Exhumed', color='blue', capsize=5)

            ax[1, 0].errorbar(x_weighted + offset, weighted_avg_maxPP_stagnant, 
                            yerr=np.array(yerr_stagnant).T, 
                            fmt='o', label='Stagnant', color='orange', capsize=5)

            # Configure ax[1, 0]
            ax[1, 0].set_xticks(x_weighted)
            ax[1, 0].set_xticklabels(wrapped_model_names, ha="center")
            ax[1, 0].set_ylabel("Weighted Average Peak P")
            ax[1, 0].set_title("Weighted Average Peak P Conditions")
            ax[1, 0].legend()

            
            ############################################################################################
            ### Scatter plot for weighted average maxPP grouped by lithology (ax[1, 1]) for exhumed particles ###
            ############################################################################################
            weighted_avg_maxPP_by_lithology = {}
            min_maxPP_by_lithology = {}
            max_maxPP_by_lithology = {}

            # Initialize dictionaries to hold maxPP values grouped by lithology
            for lithology in lithologies_exhumed_all:
                weighted_avg_maxPP_by_lithology[lithology] = []
                min_maxPP_by_lithology[lithology] = []
                max_maxPP_by_lithology[lithology] = []

            # Calculate weighted average and min/max for each lithology
            for ind_m, m in enumerate(models[test]):
                text_loc = f"{plot_loc}/{m}/txt_files"
                exhumed = pd.read_csv(f"{text_loc}/exhumed_particles.txt", sep="\t")

                # Group by lithology
                grouped = exhumed.groupby("lithology").apply(lambda x: pd.Series({
                    'weighted_mean': np.average(x['maxPP'], weights=np.ones(len(x))),  # Placeholder for actual weights if available
                    'min': np.min(x['maxPP']),
                    'max': np.max(x['maxPP'])
                })).reset_index()
                
                # Store the weighted average and min/max values for each lithology
                for index, row in grouped.iterrows():
                    lithology = row['lithology']
                    weighted_avg_maxPP_by_lithology[lithology].append(row['weighted_mean'])
                    min_maxPP_by_lithology[lithology].append(row['min'])
                    max_maxPP_by_lithology[lithology].append(row['max'])
                
                # Append NaN for missing lithologies in other models
                for lithology in lithologies_exhumed_all:
                    if lithology not in grouped['lithology'].values:
                        weighted_avg_maxPP_by_lithology[lithology].append(np.nan)
                        min_maxPP_by_lithology[lithology].append(np.nan)
                        max_maxPP_by_lithology[lithology].append(np.nan)

            # Prepare data for plotting
            offset = 0.1  # Adjust this value as needed
            for j, lithology in enumerate(lithologies_exhumed_all):
                yerr = [np.array(weighted_avg_maxPP_by_lithology[lithology]) - np.array(min_maxPP_by_lithology[lithology]), 
                        np.array(max_maxPP_by_lithology[lithology]) - np.array(weighted_avg_maxPP_by_lithology[lithology])]
                ax[1, 1].errorbar(x + j * offset, weighted_avg_maxPP_by_lithology[lithology], 
                                  yerr=yerr, 
                                  fmt='o', label=lithology, capsize=5)

            ax[1, 1].set_xticks(x)
            ax[1, 1].set_xticklabels(wrapped_model_names, ha="center")
            ax[1, 1].set_ylabel("Weighted Average MaxPP")
            ax[1, 1].set_title("Weighted Average MaxPP by Lithology")
            ax[1, 1].legend(title="Lithology")

            ############################################################################################
            ### Scatter plot for weighted average maxPP grouped by lithology (ax[1, 2]) for stagnant particles ###
            ############################################################################################
            weighted_avg_maxPP_by_lithology_stagnant = {}
            min_maxPP_by_lithology_stagnant = {}
            max_maxPP_by_lithology_stagnant = {}

            # Initialize dictionaries to hold maxPP values grouped by lithology
            for lithology in lithologies_stagnant_all:
                weighted_avg_maxPP_by_lithology_stagnant[lithology] = []
                min_maxPP_by_lithology_stagnant[lithology] = []
                max_maxPP_by_lithology_stagnant[lithology] = []

            # Calculate weighted average and min/max for each lithology
            for ind_m, m in enumerate(models[test]):
                text_loc = f"{plot_loc}/{m}/txt_files"
                stagnant = pd.read_csv(f"{text_loc}/stagnant_particles.txt", sep="\t")

                # Group by lithology
                grouped = stagnant.groupby("lithology").apply(lambda x: pd.Series({
                    'weighted_mean': np.average(x['maxPP'], weights=np.ones(len(x))),  # Placeholder for actual weights if available
                    'min': np.min(x['maxPP']),
                    'max': np.max(x['maxPP'])
                })).reset_index()
                
                # Store the weighted average and min/max values for each lithology
                for index, row in grouped.iterrows():
                    lithology = row['lithology']
                    weighted_avg_maxPP_by_lithology_stagnant[lithology].append(row['weighted_mean'])
                    min_maxPP_by_lithology_stagnant[lithology].append(row['min'])
                    max_maxPP_by_lithology_stagnant[lithology].append(row['max'])
                
                # Append NaN for missing lithologies in other models
                for lithology in lithologies_stagnant_all:
                    if lithology not in grouped['lithology'].values:
                        weighted_avg_maxPP_by_lithology_stagnant[lithology].append(np.nan)
                        min_maxPP_by_lithology_stagnant[lithology].append(np.nan)
                        max_maxPP_by_lithology_stagnant[lithology].append(np.nan)

            # Prepare data for plotting
            for j, lithology in enumerate(lithologies_stagnant_all):
                yerr = [np.array(weighted_avg_maxPP_by_lithology_stagnant[lithology]) - np.array(min_maxPP_by_lithology_stagnant[lithology]), 
                        np.array(max_maxPP_by_lithology_stagnant[lithology]) - np.array(weighted_avg_maxPP_by_lithology_stagnant[lithology])]
                ax[1, 2].errorbar(x + j * offset, weighted_avg_maxPP_by_lithology_stagnant[lithology], 
                                  yerr=yerr, 
                                  fmt='o', label=lithology, capsize=5)

            ax[1, 2].set_xticks(x)
            ax[1, 2].set_xticklabels(wrapped_model_names, ha="center")
            ax[1, 2].set_ylabel("Weighted Average MaxPP")
            ax[1, 2].set_title("Weighted Average MaxPP by Lithology")
            ax[1, 2].legend(title="Lithology")


                    

            
            





            # Plotting the figure
            plt.tight_layout()
            plt.subplots_adjust(top=0.88)  # Adjust for the main title
            plt.savefig(f"combined/avg_range/{test}_comparison_plot.png")
            plt.close()

if __name__ == "__main__":
    main()