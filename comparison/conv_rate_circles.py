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
    tests = ["friction"]
    names = ["friction_names"]

    # Set colorblind-friendly palette
    palette = sns.color_palette("colorblind")

    # Ensure the output directory exists
    output_dir = "combined/timing/by_lithology"
    os.makedirs(output_dir, exist_ok=True)

    # Loop through each test to create a grid plot
    for idx, test in enumerate(tests):

        fig, ax = plt.subplots(1, 1, figsize=(10, 7))

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

            # Create violin plots
            sns.violinplot(data=combined_data_exhumed, x='model', y='tin', hue='lithology', scale='width', palette=palette, width=0.5, legend=False, ax = ax)
            # sns.violinplot(data=combined_data_exhumed, x='model', y='tmax', hue='lithology', scale='width', palette=palette, dodge=True, ax=ax[0, 1], legend=False)
            # sns.violinplot(data=combined_data_exhumed, x='model', y='tfin', hue='lithology', scale='width', palette=palette, width=0.5, ax=ax[0, 2], legend=False)

            # sns.violinplot(data=combined_data_stagnant, x='model', y='tin', hue='lithology', scale='width', palette=palette, width=0.5, ax=ax[1, 0], legend=False)
            # sns.violinplot(data=combined_data_stagnant, x='model', y='tmax', hue='lithology', scale='width', palette=palette, dodge=True, ax=ax[1, 1], legend=False)
            # sns.violinplot(data=combined_data_stagnant, x='model', y='tfin', hue='lithology', scale='width', palette=palette, width=0.5, ax=ax[1, 2], legend=False)

            # # Adding titles and labels for each subplot
            # ax[0, 0].set_title("Exhumed Particles Timing: Tin")
            # ax[0, 0].set_xlabel("Model Name")
            # ax[0, 0].set_ylabel("Tin")
            # ax[0, 1].set_title("Exhumed Particles Timing: Tmax")
            # ax[0, 1].set_xlabel("Model Name")
            # ax[0, 1].set_ylabel("Tmax")
            # ax[0, 2].set_title("Exhumed Particles Timing: Tfin")
            # ax[0, 2].set_xlabel("Model Name")
            # ax[0, 2].set_ylabel("Tfin")

            # ax[1, 0].set_title("Stagnant Particles Timing: Tin")
            # ax[1, 0].set_xlabel("Model Name")
            # ax[1, 0].set_ylabel("Tin")
            # ax[1, 1].set_title("Stagnant Particles Timing: Tmax")
            # ax[1, 1].set_xlabel("Model Name")
            # ax[1, 1].set_ylabel("Tmax")
            # ax[1, 2].set_title("Stagnant Particles Timing: Tfin")
            # ax[1, 2].set_xlabel("Model Name")
            # ax[1, 2].set_ylabel("Tfin")

            # Wrap x-axis labels
            # for a in ax.flatten():  # Use flatten to access all axes
            labels = ax.get_xticklabels()
            wrapped_labels = [label.get_text().replace(' ', '\n') for label in labels]  # Wrap text on spaces
            ax.set_xticklabels(wrapped_labels, rotation=0, ha='right')  # Rotate labels for better visibility

            # Load conv_rate data
            conv_data = pd.read_csv(f"{plot_loc}/{models[test][0]}/txt_files/2D_v.txt", sep="\s+")
            conv_data["conv_rate"].iloc[0] = np.nan  # Avoid negative diameters
            conv_data["time"] = conv_data["time"]/1.e6

            # New x position for convergence rate circles (adjusted to the left side)
            conv_rate_x_position = -0.7  # Positioning slightly left of the violin plots

            # Scatter plot along the y-axis for each subplot
            # for i, ax_row in enumerate(ax):
            #     for j, axis in enumerate(ax_row):
                    # Determine y-axis limits for cropping
            y_limits = ax.get_ylim()

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
            ax.scatter([conv_rate_x_position] * len(cropped_conv_data), 
                            cropped_conv_data['time'], 
                            s=sizes, alpha=0.5, edgecolor='k', color = "grey")  # Red circles for convergence rate
                    

            # Save the plot
            plt.tight_layout()
            plt.show()
            # plt.savefig(f"{output_dir}/exhumed_timing_combined_{test}.png")
            # plt.close()  # Close the figure after saving

if __name__ == "__main__":
    main()
