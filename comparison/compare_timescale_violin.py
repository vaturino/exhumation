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
    lithology_colors = {
        'sed': palette[0],
        'oc': palette[1],
        'ecl': palette[2],
        'serp': palette[3],
        # Add more lithologies as needed
    }

    # Ensure the output directory exists
    output_dir = "combined/timing/by_lithology"
    os.makedirs(output_dir, exist_ok=True)
    
    # Create subdirectories for exhumed and stagnant data
    exhumed_dir = os.path.join(output_dir, "exhumed")
    stagnant_dir = os.path.join(output_dir, "stagnant")
    os.makedirs(exhumed_dir, exist_ok=True)
    os.makedirs(stagnant_dir, exist_ok=True)

    # Loop through each test to create a grid plot
    for idx, test in enumerate(tests):
        # Initialize separate plots for exhumed and stagnant particles
        fig_exhumed, ax_exhumed = plt.subplots(1, 3, figsize=(17, 5))
        fig_stagnant, ax_stagnant = plt.subplots(1, 3, figsize=(17, 5))

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

            # print(combined_data_stagnant[(combined_data_stagnant['model'] == model_names[0]) & (combined_data_stagnant['lithology'] == 'ecl')])  
            # exit()

            # Create violin plots for exhumed particles
            sns.violinplot(data=combined_data_exhumed, x='model', y='tin', hue='lithology', scale='width', palette=lithology_colors, width=0.5, ax=ax_exhumed[0], legend=False)
            sns.violinplot(data=combined_data_exhumed, x='model', y='tmax', hue='lithology', scale='width', palette=lithology_colors, dodge=True, ax=ax_exhumed[1], legend=False)
            sns.violinplot(data=combined_data_exhumed, x='model', y='tfin', hue='lithology', scale='width', palette=lithology_colors, width=0.5, ax=ax_exhumed[2], legend=False)

            # Create violin plots for stagnant particles
            sns.violinplot(data=combined_data_stagnant, x='model', y='tin', hue='lithology', scale='width', palette=lithology_colors, width=0.5, ax=ax_stagnant[0], legend=False)
            sns.violinplot(data=combined_data_stagnant, x='model', y='tmax', hue='lithology', scale='width', palette=lithology_colors, dodge=True, ax=ax_stagnant[1], legend=False)
            sns.violinplot(data=combined_data_stagnant, x='model', y='tfin', hue='lithology', scale='width', palette=lithology_colors, width=0.5, ax=ax_stagnant[2], legend=False)

            # Adding titles and labels for exhumed particles
            ax_exhumed[0].set_title("Exhumed Particles Timing: Tin")
            ax_exhumed[0].set_xlabel("Model Name")
            ax_exhumed[0].set_ylabel("Tin")
            ax_exhumed[1].set_title("Exhumed Particles Timing: Tmax")
            ax_exhumed[1].set_xlabel("Model Name")
            ax_exhumed[1].set_ylabel("Tmax")
            ax_exhumed[2].set_title("Exhumed Particles Timing: Tfin")
            ax_exhumed[2].set_xlabel("Model Name")
            ax_exhumed[2].set_ylabel("Tfin")

            # Adding titles and labels for stagnant particles
            ax_stagnant[0].set_title("Stagnant Particles Timing: Tin")
            ax_stagnant[0].set_xlabel("Model Name")
            ax_stagnant[0].set_ylabel("Tin")
            ax_stagnant[1].set_title("Stagnant Particles Timing: Tmax")
            ax_stagnant[1].set_xlabel("Model Name")
            ax_stagnant[1].set_ylabel("Tmax")
            ax_stagnant[2].set_title("Stagnant Particles Timing: Tfin")
            ax_stagnant[2].set_xlabel("Model Name")
            ax_stagnant[2].set_ylabel("Tfin")

            # Wrap x-axis labels
            for a in ax_exhumed:  # Use exhumed axes
                labels = a.get_xticklabels()
                wrapped_labels = [label.get_text().replace(' ', '\n') for label in labels]  # Wrap text on spaces
                a.set_xticklabels(wrapped_labels, rotation=0, ha='right')  # Rotate labels for better visibility

            for a in ax_stagnant:  # Use stagnant axes
                labels = a.get_xticklabels()
                wrapped_labels = [label.get_text().replace(' ', '\n') for label in labels]  # Wrap text on spaces
                a.set_xticklabels(wrapped_labels, rotation=0, ha='right')  # Rotate labels for better visibility

            # Load conv_rate data
            conv_data = pd.read_csv(f"{plot_loc}/{models[test][0]}/txt_files/2D_v.txt", sep="\s+")
            conv_data["conv_rate"].iloc[0] = np.nan  # Avoid negative diameters
            conv_data["time"] = conv_data["time"]/1.e6

            # New x position for convergence rate circles (adjusted to the left side)
            conv_rate_x_position = -0.7  # Positioning slightly left of the violin plots

            # Scatter plot along the y-axis for each subplot
            for i, ax_row in enumerate([ax_exhumed, ax_stagnant]):
                for j, axis in enumerate(ax_row):
                    # Determine y-axis limits for cropping
                    y_limits = axis.get_ylim()

                    # Crop conv_data to only include time values within the y-axis limits
                    cropped_conv_data = conv_data[(conv_data['time'] >= y_limits[0]) & (conv_data['time'] <= y_limits[1])]

                    # Normalize conv_rate for better visualization of circle sizes
                    min_conv_rate = cropped_conv_data['conv_rate'].min()
                    adjusted_conv_rate = cropped_conv_data['conv_rate'] - min_conv_rate
                    adjusted_conv_rate[adjusted_conv_rate < 0] = 0  # Avoid negative diameters

                    # Set circle sizes (diameter) proportional to conv_rate
                    size_factor = 100  # Adjust this factor as needed
                    sizes = adjusted_conv_rate * size_factor

                    # Scatter plot for convergence rate circles on the left side
                    axis.scatter([conv_rate_x_position] * len(cropped_conv_data), 
                                 cropped_conv_data['time'], 
                                 s=sizes, alpha=0.5, edgecolor='k', color="grey")  # Red circles for convergence rate

            # Save the plots
            plt.tight_layout()
            fig_exhumed.savefig(f"{exhumed_dir}/exhumed_timing_combined_{test}.png")
            plt.close(fig_exhumed)  # Close the exhumed figure after saving
            plt.tight_layout()
            fig_stagnant.savefig(f"{stagnant_dir}/stagnant_timing_combined_{test}.png")
            plt.close(fig_stagnant)  # Close the stagnant figure after saving

if __name__ == "__main__":
    main()
