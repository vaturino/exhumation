#!/usr/bin/env python3

import pandas as pd
import seaborn as sns
import json
import os
import matplotlib.pyplot as plt
import joypy  # Ensure you have the joypy library installed
import numpy as np
from PIL import Image
import matplotlib.gridspec as gridspec

def main():
    # Load json file
    with open('models.json') as f:
        models = json.load(f)

    plot_loc = "/home/vturino/PhD/projects/exhumation/plots/single_models"
    tests = ["velocity", "viscosity", "friction", "serpentinization"]
    names = ["velocity_names", "viscosity_names", "friction_names", "serpentinization_names"]

    # Set colorblind-friendly palette
    palette = sns.color_palette("colorblind")

    # Ensure the output directory exists
    output_dir = "combined/timing/by_lithology"
    os.makedirs(output_dir, exist_ok=True)

    # Loop through each test to create a grid plot
    for idx, test in enumerate(tests):
        if test in models:
            model_names = models[names[tests.index(test)]]

            # Prepare data for Ridge Plots
            ridge_data_tin_exhumed = []
            ridge_data_tmax_exhumed = []
            ridge_data_tfin_exhumed = []
            labels_tin_exhumed = []
            labels_tmax_exhumed = []
            labels_tfin_exhumed = []
            
            ridge_data_tin_stagnant = []
            ridge_data_tmax_stagnant = []
            ridge_data_tfin_stagnant = []
            labels_tin_stagnant = []
            labels_tmax_stagnant = []
            labels_tfin_stagnant = []
            
            lithology_colors = {}

            # Loop through each model to collect data
            for ind_m, m in enumerate(models[test]):
                text_loc = f"{plot_loc}/{m}/txt_files"
                
                try:
                    exhumed = pd.read_csv(f"{text_loc}/timing_exhumed_particles.txt", sep="\s+")
                    stagnant = pd.read_csv(f"{text_loc}/timing_stagnant_particles.txt", sep="\s+")
                except FileNotFoundError:
                    print(f"Files not found for model {m}. Skipping.")
                    continue

                # Assign colors to lithologies if not already done
                for lith in exhumed["lithology"].unique():
                    if lith not in lithology_colors:
                        lithology_colors[lith] = palette[len(lithology_colors) % len(palette)]
                
                for lith in stagnant["lithology"].unique():
                    if lith not in lithology_colors:
                        lithology_colors[lith] = palette[len(lithology_colors) % len(palette)]

                # Collecting data for Ridge Plot (exhumed)
                for lith in exhumed["lithology"].unique():
                    lith_data_tin = exhumed[exhumed["lithology"] == lith]["tin"].values
                    lith_data_tmax = exhumed[exhumed["lithology"] == lith]["tmax"].values
                    lith_data_tfin = exhumed[exhumed["lithology"] == lith]["tfin"].values
                    
                    if lith_data_tin.size > 0:
                        ridge_data_tin_exhumed.append(lith_data_tin)
                        labels_tin_exhumed.append(f"{model_names[ind_m]} - {lith}")  # Label for tin
                    if lith_data_tmax.size > 0:
                        ridge_data_tmax_exhumed.append(lith_data_tmax)
                        labels_tmax_exhumed.append(f"{model_names[ind_m]} - {lith}")  # Label for tmax
                    if lith_data_tfin.size > 0:
                        ridge_data_tfin_exhumed.append(lith_data_tfin)
                        labels_tfin_exhumed.append(f"{model_names[ind_m]} - {lith}")  # Label for tfin

                # Collecting data for Ridge Plot (stagnant)
                for lith in stagnant["lithology"].unique():
                    lith_data_tin = stagnant[stagnant["lithology"] == lith]["tin"].values
                    lith_data_tmax = stagnant[stagnant["lithology"] == lith]["tmax"].values
                    lith_data_tfin = stagnant[stagnant["lithology"] == lith]["tfin"].values
                    
                    if lith_data_tin.size > 0:
                        ridge_data_tin_stagnant.append(lith_data_tin)
                        labels_tin_stagnant.append(f"{model_names[ind_m]} - {lith}")  # Label for tin
                    if lith_data_tmax.size > 0:
                        ridge_data_tmax_stagnant.append(lith_data_tmax)
                        labels_tmax_stagnant.append(f"{model_names[ind_m]} - {lith}")  # Label for tmax
                    if lith_data_tfin.size > 0:
                        ridge_data_tfin_stagnant.append(lith_data_tfin)
                        labels_tfin_stagnant.append(f"{model_names[ind_m]} - {lith}")  # Label for tfin

            # Pad data to ensure all arrays have the same length
            def pad_data(data):
                if not data:
                    return []
                max_len = max(len(arr) for arr in data)
                return [np.pad(arr, (0, max_len - len(arr)), constant_values=np.nan) for arr in data]

            ridge_data_tin_exhumed = pad_data(ridge_data_tin_exhumed)
            ridge_data_tmax_exhumed = pad_data(ridge_data_tmax_exhumed)
            ridge_data_tfin_exhumed = pad_data(ridge_data_tfin_exhumed)

            ridge_data_tin_stagnant = pad_data(ridge_data_tin_stagnant)
            ridge_data_tmax_stagnant = pad_data(ridge_data_tmax_stagnant)
            ridge_data_tfin_stagnant = pad_data(ridge_data_tfin_stagnant)

            # Determine height ratios based on unique lithologies
            unique_exhumed_lithologies = exhumed["lithology"].unique()
            unique_stagnant_lithologies = stagnant["lithology"].unique()

            # Adjust height ratios: taller for more stagnant lithologies
            height_ratios = [1, len(unique_stagnant_lithologies)]  # First row normal height, second row based on stagnant lithologies

            # Create a final figure for this test with gridspec for dynamic row heights
            fig = plt.figure(figsize=(15, 8))
            gs = gridspec.GridSpec(2, 3)  # 2 rows, 3 columns with custom heights

            # Create and save the joy plots for exhumed particles
            def create_joy_plot(data, labels, title, output_file, y_scale=1):
                joy_fig, joy_ax = plt.subplots(figsize=(8, 6))
                joypy.joyplot(data=data, overlap=0.5, labels=labels,
                            color=[lithology_colors[label.split(" - ")[1]] for label in labels], ax=joy_ax)

                # Set x-axis limits to be consistent across all plots
                joy_ax.set_xlim(0, 50)  # Adjust this based on your actual data range
                joy_ax.set_title(title, fontsize=20)  # Title font size
                joy_ax.set_xlabel("Subduction time (Ma)", fontsize=20)  # X-axis label font size
                joy_ax.set_ylabel("Density", fontsize=20)  # Y-axis label font size
                joy_ax.tick_params(axis='both', which='major', labelsize=20)  # Tick label size

                # Apply vertical scaling to the y-axis for stagnant plots
                if y_scale != 1:
                    y_ticks = joy_ax.get_yticks()
                    joy_ax.set_ylim(0, y_ticks[-1] * y_scale)  # Increase the upper limit based on y_scale

                plt.savefig(output_file)
                plt.close()  # Close the figure after saving to free up memory


            # Create joy plots for exhumed particles
            joy_output_file_tin_exhumed = f"{output_dir}/single_exhumed/ridge_plot_tin_{test}.png"
            create_joy_plot(ridge_data_tin_exhumed, labels_tin_exhumed, f"Ridge Plot of Exhumed Particles by Lithology and Model (tin) - {test}", joy_output_file_tin_exhumed)
            ax_tin_exhumed = fig.add_subplot(gs[0, 0])
            ax_tin_exhumed.imshow(Image.open(joy_output_file_tin_exhumed))
            ax_tin_exhumed.axis('off')
            ax_tin_exhumed.set_title(f"{test} - tin", fontsize=14)

            joy_output_file_tmax_exhumed = f"{output_dir}/single_exhumed/ridge_plot_tmax_{test}.png"
            create_joy_plot(ridge_data_tmax_exhumed, labels_tmax_exhumed, f"Ridge Plot of Exhumed Particles by Lithology and Model (tmax) - {test}", joy_output_file_tmax_exhumed)
            ax_tmax_exhumed = fig.add_subplot(gs[0, 1])
            ax_tmax_exhumed.imshow(Image.open(joy_output_file_tmax_exhumed))
            ax_tmax_exhumed.axis('off')
            ax_tmax_exhumed.set_title(f"{test} - tmax", fontsize=14)

            joy_output_file_tfin_exhumed = f"{output_dir}/single_exhumed/ridge_plot_tfin_{test}.png"
            create_joy_plot(ridge_data_tfin_exhumed, labels_tfin_exhumed, f"Ridge Plot of Exhumed Particles by Lithology and Model (tfin) - {test}", joy_output_file_tfin_exhumed)
            ax_tfin_exhumed = fig.add_subplot(gs[0, 2])
            ax_tfin_exhumed.imshow(Image.open(joy_output_file_tfin_exhumed))
            ax_tfin_exhumed.axis('off')
            ax_tfin_exhumed.set_title(f"{test} - tfin", fontsize=14)

            # Create joy plots for stagnant particles
            joy_output_file_tin_stagnant = f"{output_dir}/single_stagnant/ridge_plot_tin_{test}.png"
            create_joy_plot(ridge_data_tin_stagnant, labels_tin_stagnant, f"Ridge Plot of Stagnant Particles by Lithology and Model (tin) - {test}", joy_output_file_tin_stagnant)
            ax_tin_stagnant = fig.add_subplot(gs[1, 0])
            ax_tin_stagnant.imshow(Image.open(joy_output_file_tin_stagnant))
            ax_tin_stagnant.axis('off')
            ax_tin_stagnant.set_title(f"{test} - tin (stagnant)", fontsize=14)

            joy_output_file_tmax_stagnant = f"{output_dir}/single_stagnant/ridge_plot_tmax_{test}.png"
            create_joy_plot(ridge_data_tmax_stagnant, labels_tmax_stagnant, f"Ridge Plot of Stagnant Particles by Lithology and Model (tmax) - {test}", joy_output_file_tmax_stagnant)
            ax_tmax_stagnant = fig.add_subplot(gs[1, 1])
            ax_tmax_stagnant.imshow(Image.open(joy_output_file_tmax_stagnant))
            ax_tmax_stagnant.axis('off')
            ax_tmax_stagnant.set_title(f"{test} - tmax (stagnant)", fontsize=14)

            joy_output_file_tfin_stagnant = f"{output_dir}/single_stagnant/ridge_plot_tfin_{test}.png"
            create_joy_plot(ridge_data_tfin_stagnant, labels_tfin_stagnant, f"Ridge Plot of Stagnant Particles by Lithology and Model (tfin) - {test}", joy_output_file_tfin_stagnant)
            ax_tfin_stagnant = fig.add_subplot(gs[1, 2])
            ax_tfin_stagnant.imshow(Image.open(joy_output_file_tfin_stagnant))
            ax_tfin_stagnant.axis('off')
            ax_tfin_stagnant.set_title(f"{test} - tfin (stagnant)", fontsize=14)

            # Adjust layout for better spacing and save the final figure for this test
            plt.tight_layout()
            final_output_file = f"{output_dir}/ridge_plot_grid_{test}.png"
            plt.savefig(final_output_file)
            plt.close()  # Close the figure after saving to free up memory

if __name__ == "__main__":
    main()
