#!/usr/bin/env python3

import pandas as pd
import seaborn as sns
import json
import os
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches


def draw_rectangles(ax, data, x, y, hue, colors, alpha, size):
    """Draw rectangular markers for stripplot."""
    for i, (name, group) in enumerate(data.groupby(hue)):
        x_pos = np.full(group[x].shape, i)  # x position for each group
        # Draw rectangles with height < width
        for j in range(group.shape[0]):
            rect = mpatches.Rectangle((x_pos[j] - 0.25, group[y].iloc[j]), 0.5, 0.2, 
                                      color=colors[name], alpha=alpha)
            ax.add_patch(rect)

def main():
    # Load json file
    with open('models.json') as f:
        models = json.load(f)

    plot_loc = "/home/vturino/PhD/projects/exhumation/plots/single_models"
    tests = ["velocity", "viscosity", "friction", "serpentinization"]
    names = ["velocity_names", "viscosity_names", "friction_names", "serpentinization_names"]

    # Set colorblind-friendly palette
    palette = sns.color_palette("colorblind")

    # Define color palettes
    colors_tin = {
        "sed": "midnightblue",
        "oc": "#733A11",
        "ecl": "darkgreen",
        "serp": "#3b0000"
    }

    colors_tmax = {
        "sed": "mediumblue",
        "oc": "#B06D1A",
        "ecl": "forestgreen",
        "serp": "brown"
    }

    colors_tfin = {
        "sed": "cornflowerblue",
        "oc": "#E3B64F",
        "ecl": "olivedrab",
        "serp": "lightsalmon"
    }

    # Ensure the output directory exists
    output_dir = "combined/timing/strip_plots"
    os.makedirs(output_dir, exist_ok=True)

    alpha = 0.15
    size = 15

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

            # Call the custom rectangle drawing function for exhumed particles
            draw_rectangles(ax[0], combined_data_exhumed, 'model', 'tin', 'lithology', colors_tin, alpha, size)
            draw_rectangles(ax[0], combined_data_exhumed, 'model', 'tmax', 'lithology', colors_tmax, alpha, size)
            draw_rectangles(ax[0], combined_data_exhumed, 'model', 'tfin', 'lithology', colors_tfin, alpha, size)

            # Set the x-axis labels and title
            ax[0].set_xlabel('Model')
            ax[0].set_ylabel('Timing Values')
            ax[0].set_title(f'Exhumed Particles Timing - {test}')

            # (Optional) Create a similar plot for stagnant particles if needed
            # Uncomment below to draw for stagnant particles
            # draw_rectangles(ax[1], combined_data_stagnant, 'model', 'tin', 'lithology', colors_tin, alpha, size)
            # draw_rectangles(ax[1], combined_data_stagnant, 'model', 'tmax', 'lithology', colors_tmax, alpha, size)
            # draw_rectangles(ax[1], combined_data_stagnant, 'model', 'tfin', 'lithology', colors_tfin, alpha, size)
            # ax[1].set_xlabel('Model')
            # ax[1].set_ylabel('Timing Values')
            # ax[1].set_title(f'Stagnant Particles Timing - {test}')

            plt.show()
            exit()
if __name__ == "__main__":
    main()
