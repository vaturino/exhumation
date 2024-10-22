#!/usr/bin/env python3

import pandas as pd
import seaborn as sns
import json
import os
import matplotlib.pyplot as plt
import numpy as np

def main():
    # Load json file
    with open('models.json') as f:
        models = json.load(f)

    plot_loc = "/home/vturino/PhD/projects/exhumation/plots/single_models"
    tests = ["velocity", "viscosity", "friction", "serpentinization"]
    names = ["velocity_names", "viscosity_names", "friction_names", "serpentinization_names"]

    # Set colorblind-friendly palette
    palette = sns.color_palette("colorblind")

    # Create a mapping for lithology colors
     # Define color palettes
    colors_tin = {
        "sed": "midnightblue",
        "oc": "saddlebrown",
        "ecl": "darkgreen",
        "serp": "maroon"
    }

    colors_tmax = {
        "sed": "mediumblue",
        "oc": "darkorange",
        "ecl": "forestgreen",
        "serp": "firebrick"
    }

    colors_tfin = {
        "sed": "cornflowerblue",
        "oc": "goldenrod",
        "ecl": "olivedrab",
        "serp": "indianred"
    }


    # Ensure the output directory exists
    output_dir = "combined/timing/boxen_plots"
    os.makedirs(output_dir, exist_ok=True)

    k = 0.1

    # Loop through each test to create a grid plot
    for idx, test in enumerate(tests):
        fig, ax = plt.subplots(1, 2, figsize=(17, 8))

        if test in models:
            model_names = models[names[tests.index(test)]]
            all_data_exhumed = []  # List to store all dataframes for exhumed particles
            all_data_stagnant = []  # List to store all dataframes for stagnant particles

            for ind_m, m in enumerate(models[test]):
                text_loc = f"{plot_loc}/{m}/txt_files"
                exhumed = pd.read_csv(f"{text_loc}/timing_exhumed_particles.txt", sep="\s+")
                stagnant = pd.read_csv(f"{text_loc}/timing_stagnant_particles.txt", sep="\s+")

                # Add a new column for the model name from model_names
                exhumed['model'] = model_names[ind_m]
                all_data_exhumed.append(exhumed)
                stagnant['model'] = model_names[ind_m]
                all_data_stagnant.append(stagnant)

            # Combine all data into a single DataFrame
            combined_data_exhumed = pd.concat(all_data_exhumed)
            combined_data_stagnant = pd.concat(all_data_stagnant)

            alpha = 0.5
            size = 6

            sns.boxenplot(data=combined_data_exhumed,
               x='model', 
               y='tin',
               hue='lithology', 
               palette=colors_tin,
               ax=ax[0],
               dodge=True,
               linewidth=1,
               saturation=.9,
               gap=0.2,
               width=0.5,
               showfliers=True,
               width_method='linear',
               line_kws=dict(linewidth=1.5, color="#cde"))
        
            sns.boxenplot(data=combined_data_exhumed,
                x='model', 
                y='tmax',
                hue='lithology', 
                palette=colors_tmax,
                ax=ax[0],
                dodge=True,
                linewidth=1,
                saturation=.9,
                gap=0.2,
                width=0.5,
                showfliers=True,
                width_method='linear',
                line_kws=dict(linewidth=1.5, color="#cde"))
            
            sns.boxenplot(data=combined_data_exhumed,
                x='model', 
                y='tfin',
                hue='lithology', 
                palette=colors_tfin,
                ax=ax[0],
                dodge=True,
                linewidth=1,
                saturation=.9,
                gap=0.2,
                width=0.5,
                showfliers=True,
                width_method='linear',
                line_kws=dict(linewidth=1.5, color="#cde"))
            
            sns.boxenplot(data=combined_data_stagnant,
                x='model', 
                y='tin',
                hue='lithology', 
                palette=colors_tin,
                ax=ax[1],
                dodge=True,
                linewidth=1,
                saturation=.9,
                gap=0.2,
                width=0.5,
                showfliers=True,
                width_method='linear',
                line_kws=dict(linewidth=1.5, color="#cde"))
            
            sns.boxenplot(data=combined_data_stagnant,
                x='model', 
                y='tmax',
                hue='lithology', 
                palette=colors_tmax,
                ax=ax[1],
                dodge=True,
                linewidth=1,
                saturation=.9,
                gap=0.2,
                width=0.5,
                showfliers=True,
                width_method='linear',
                line_kws=dict(linewidth=1.5, color="#cde"))
            
            sns.boxenplot(data=combined_data_stagnant,
                x='model', 
                y='tfin',
                hue='lithology', 
                palette=colors_tfin,
                ax=ax[1],
                dodge=True,
                linewidth=1,
                saturation=.9,
                gap=0.2,
                width=0.5,
                showfliers=True,
                width_method='linear',
                line_kws=dict(linewidth=1.5, color="#cde"))

            # Load conv_rate data
            conv_data = pd.read_csv(f"{plot_loc}/{models[test][0]}/txt_files/2D_v.txt", sep="\s+")
            conv_data["conv_rate"].iloc[0] = np.nan  # Avoid negative diameters
            conv_data["time"] = conv_data["time"]/1.e6

            # New x position for convergence rate circles (adjusted to the left side)
            conv_rate_x_position = -0.5  # Positioning slightly left of the violin plots

            # Scatter plot along the y-axis for each subplot
            for j, axis in enumerate(ax):
                # Determine y-axis limits for cropping
                y_limits = axis.get_ylim()

                # Set circle sizes (diameter) proportional to conv_rate
                size_factor = 100  # Adjust this factor as needed
                sizes = conv_data["conv_rate"] * size_factor

                # Scatter plot for convergence rate circles on the left side
                axis.scatter([conv_rate_x_position] * len(conv_data["conv_rate"]), 
                                conv_data['time'], 
                                s=sizes, alpha=0.5, edgecolor='k', color="grey")  # Grey circles for convergence rate
                
            
            # Adding titles and labels for each subplot
            ax[0].set_title("Exhumed Particles Timing")
            ax[0].set_xlabel("Model Name")
            ax[0].set_ylabel("Time (Myr)")
            ax[0].set_ylim(0, 50)

            ax[1].set_title("Stagnant Particles Timing")
            ax[1].set_xlabel("Model Name")
            # ax[1].set_ylabel("Time (Myr)")
            ax[1].set_ylabel("")
            ax[1].set_ylim(0, 50)


            # # Create the legend for ax[1]
            # legend1 = ax[1].legend(loc='best', bbox_to_anchor=(1, 1), handletextpad=1.5, labelspacing=1)
            # legend1.set_title("Lithology")

            # Get handles and labels from ax[1]
            handles, labels = ax[1].get_legend_handles_labels()

            # Add the legend to ax[0]
            ax[0].legend(handles, labels, loc='best', bbox_to_anchor=(1, 1), handletextpad=1.5, labelspacing=1, title="Lithology")
            
            # Save the plot
            plt.tight_layout()
            format = ["png", "eps"]
            for form in format:
                plt.savefig(f"{output_dir}/{form}/timing_combined_{test}.{form}", dpi = 500, format = form)
            # plt.savefig(f"{output_dir}/exhumed_timing_combined_{test}.eps", dpi = 500, format = "eps")
            plt.close()  # Close the figure after saving

if __name__ == "__main__":
    main()
