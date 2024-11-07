#!/usr/bin/env python3

import pandas as pd
import seaborn as sns
import json
import os
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches


def main():
    # Load json file
    with open('models.json') as f:
        models = json.load(f)

    plot_loc = "/home/vturino/PhD/projects/exhumation/plots/single_models"
    tests = ["velocity"] #, "viscosity", "friction", "serpentinization"]
    names = ["velocity_names"] #, "viscosity_names", "friction_names", "serpentinization_names"]

    # Create a mapping for lithology colors
     # Define color palettes
    colors_tin = {
        "sed": "midnightblue",
        "oc": "#733A11",
        "ecl": "#003300",
        "serp": "#3b0000"
    }

    colors_tmax = {
        "sed": "mediumblue",
        "oc": "#B06D1A",
        "ecl": "#45701C",
        "serp": "brown"
    }

    colors_tfin = {
        "sed": "cornflowerblue",
        "oc": "#E3B64F",
        "ecl": "#A0C93D",
        "serp": "lightsalmon"
    }


    # Ensure the output directory exists
    output_dir = "combined/timing/strip_plots"
    os.makedirs(output_dir, exist_ok=True)

    # Loop through each test to create a grid plot
    for idx, test in enumerate(tests):
        fig, ax = plt.subplots(1, 2, figsize=(20, 10), gridspec_kw={'width_ratios': [1, 1.5]})

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
                exhumed[f"{test}_values"] = models[f"{test}_values"][ind_m]
                stagnant[f"{test}_values"] = models[f"{test}_values"][ind_m]

            # Combine all data into a single DataFrame
            combined_data_exhumed = pd.concat(all_data_exhumed)
            combined_data_stagnant = pd.concat(all_data_stagnant)


            combined_data_exhumed = combined_data_exhumed[combined_data_exhumed[f'{test}_values'] != 10]
            combined_data_stagnant = combined_data_stagnant[combined_data_stagnant[f'{test}_values'] != 10]

            


            # CREATE SCATTER PLOT

            ########### EXHUMED PARTICLES ###########
            # Group by 'tmax' and 'lithology' to get the count of each combination
            tmax_lith_counts_exhumed = combined_data_exhumed.groupby(['tmax', 'lithology', f'{test}_values']).size().unstack(fill_value=0)
            # Convert the DataFrame to a long format
            long_df_exhumed = tmax_lith_counts_exhumed.reset_index().melt(id_vars=['tmax', 'lithology'], var_name=f'{test}_values', value_name='count')
            exhumed_scale = 1
            for l, lith in enumerate(long_df_exhumed['lithology'].unique()):
                exhumed_scale = long_df_exhumed[long_df_exhumed['lithology'] == lith]['count'].max()
                # exhumed_scale the lith count by the maximum count across all lithologies
                long_df_exhumed.loc[long_df_exhumed['lithology'] == lith, 'count'] = long_df_exhumed[long_df_exhumed['lithology'] == lith]['count'] / exhumed_scale * 100

            dodge_offsets_exhumed = {hue: i * 3 for i, hue in enumerate(long_df_exhumed['lithology'].unique())}
            long_df_exhumed['x_dodge'] = long_df_exhumed.apply(
                lambda row: row[f'{test}_values'] + dodge_offsets_exhumed[row['lithology']], axis=1
            )


            ########### STAGNANT PARTICLES ###########
            # Group by 'tmax' and 'lithology' to get the count of each combination
            tmax_lith_counts_stagnant = combined_data_stagnant.groupby(['tmax', 'lithology', f'{test}_values']).size().unstack(fill_value=0)
            # Convert the DataFrame to a long format
            long_df_stagnant = tmax_lith_counts_stagnant.reset_index().melt(id_vars=['tmax', 'lithology'], var_name=f'{test}_values', value_name='count')
            stagnant_scale = 1
            for l, lith in enumerate(long_df_stagnant['lithology'].unique()):
                stagnant_scale = long_df_stagnant[long_df_stagnant['lithology'] == lith]['count'].max()
                # exhumed_scale the lith count by the maximum count across all lithologies
                long_df_stagnant.loc[long_df_stagnant['lithology'] == lith, 'count'] = long_df_stagnant[long_df_stagnant['lithology'] == lith]['count'] / stagnant_scale * 100

            dodge_offsets_stagnant = {hue: i * 5 for i, hue in enumerate(long_df_stagnant['lithology'].unique())}
            long_df_stagnant['x_dodge'] = long_df_stagnant.apply(
                lambda row: row[f'{test}_values'] + dodge_offsets_stagnant[row['lithology']], axis=1
            )



            ################### CALCULATE FOR TFIN ###################
            # Exhumed particles
            # Group by 'tfin' and 'lithology' to get the count of each combination
            tfin_lith_counts_exhumed = combined_data_exhumed.groupby(['tfin', 'lithology', f'{test}_values']).size().unstack(fill_value=0)
            # Convert the DataFrame to a long format
            long_df_exhumed_tfin = tfin_lith_counts_exhumed.reset_index().melt(id_vars=['tfin', 'lithology'], var_name=f'{test}_values', value_name='count')
            exhumed_scale_tfin = 1
            for l, lith in enumerate(long_df_exhumed_tfin['lithology'].unique()):
                exhumed_scale_tfin = long_df_exhumed_tfin[long_df_exhumed_tfin['lithology'] == lith]['count'].max()
                # exhumed_scale the lith count by the maximum count across all lithologies
                long_df_exhumed_tfin.loc[long_df_exhumed_tfin['lithology'] == lith, 'count'] = long_df_exhumed_tfin[long_df_exhumed_tfin['lithology'] == lith]['count'] / exhumed_scale_tfin * 100
            
            dodge_offsets_exhumed_tfin = {hue: i * 3 for i, hue in enumerate(long_df_exhumed_tfin['lithology'].unique())}
            long_df_exhumed_tfin['x_dodge'] = long_df_exhumed_tfin.apply(
                lambda row: row[f'{test}_values'] + dodge_offsets_exhumed_tfin[row['lithology']], axis=1
            )

            # Stagnant particles
            # Group by 'tfin' and 'lithology' to get the count of each combination
            tfin_lith_counts_stagnant = combined_data_stagnant.groupby(['tfin', 'lithology', f'{test}_values']).size().unstack(fill_value=0)
            # Convert the DataFrame to a long format
            long_df_stagnant_tfin = tfin_lith_counts_stagnant.reset_index().melt(id_vars=['tfin', 'lithology'], var_name=f'{test}_values', value_name='count')
            stagnant_scale_tfin = 1
            for l, lith in enumerate(long_df_stagnant_tfin['lithology'].unique()):
                stagnant_scale_tfin = long_df_stagnant_tfin[long_df_stagnant_tfin['lithology'] == lith]['count'].max()
                # exhumed_scale the lith count by the maximum count across all lithologies
                long_df_stagnant_tfin.loc[long_df_stagnant_tfin['lithology'] == lith, 'count'] = long_df_stagnant_tfin[long_df_stagnant_tfin['lithology'] == lith]['count'] / stagnant_scale_tfin * 100

            dodge_offsets_stagnant_tfin = {hue: i * 5 for i, hue in enumerate(long_df_stagnant_tfin['lithology'].unique())}
            long_df_stagnant_tfin['x_dodge'] = long_df_stagnant_tfin.apply(
                lambda row: row[f'{test}_values'] + dodge_offsets_stagnant_tfin[row['lithology']], axis=1
            )



            ############# SCATTER PLOT #############
           
            alpha = 1
            # marker = '$\u25AC$'
            marker = 'o'

            rect = mpatches.Rectangle((-20, 0), 200, 10, linewidth=1, edgecolor='grey', facecolor='powderblue', alpha=0.3, zorder = 1)
            ax[0].add_patch(rect)

                        
            # Now plot with the new x_dodge column
            sns.scatterplot(
                data=long_df_exhumed,
                x='x_dodge', y='tmax',
                hue='lithology', palette=colors_tmax,
                size='count', sizes=(50, 500),  # Adjust (min, max) sizes as needed
                alpha=alpha, ax=ax[0],
                linewidth=0,
                legend=False, marker = marker
            )

            ax[0].set_ylim(0, 50)
            ax[0].set_xlim(- 0.1*long_df_exhumed["x_dodge"].max(), long_df_exhumed["x_dodge"].max()+ 0.2*long_df_exhumed["x_dodge"].max())

            # Set x-axis ticks just for values of `{test}_values`
            unique_test_values = combined_data_exhumed[f'{test}_values'].unique()
            ax[0].set_xticks(unique_test_values)
            ax[0].set_xticklabels(combined_data_exhumed['model'].unique())
            


            # Now plot with the new x_dodge column
            sns.scatterplot(
                data=long_df_stagnant,
                x='x_dodge', y='tmax',
                hue='lithology', palette=colors_tmax,
                size='count', sizes=(50, 500),  # Adjust (min, max) sizes as needed
                alpha=alpha, ax=ax[1],
                linewidth=0,
                legend=False, marker = marker
            )

           
            rect2 = mpatches.Rectangle((-20, 0), 200, 10, linewidth=1, edgecolor='grey', facecolor='powderblue', alpha=0.3, zorder = 1)
            ax[1].add_patch(rect2)
            ax[1].set_xlim(- 0.1*long_df_exhumed["x_dodge"].max(), long_df_exhumed["x_dodge"].max()+ 0.2*long_df_exhumed["x_dodge"].max())
            ax[1].set_ylim(0, 50)
            ax[1].set_xticks(unique_test_values)
            ax[1].set_xticklabels(combined_data_exhumed['model'].unique())


            # # # plot tfin: dodge this by 5
            long_df_exhumed_tfin["x_fin"] = long_df_exhumed_tfin['x_dodge'] + 8
            long_df_stagnant_tfin["x_fin"] = long_df_stagnant_tfin['x_dodge'] + 8

            sns.scatterplot(
                data=long_df_exhumed_tfin,
                x='x_fin', y='tfin',
                hue='lithology', palette=colors_tfin,
                size='count', sizes=(50, 500),  # Adjust (min, max) sizes as needed
                alpha=alpha, ax=ax[0],
                linewidth=0,
                legend=False, marker = marker
            )

            sns.scatterplot(
                data=long_df_stagnant_tfin,
                x='x_fin', y='tfin',
                hue='lithology', palette=colors_tfin,
                size='count', sizes=(50, 500),  # Adjust (min, max) sizes as needed
                alpha=alpha, ax=ax[1],
                linewidth=0,
                legend=False, marker = marker
            )

             # Load conv_rate data
            conv_data = pd.read_csv(f"{plot_loc}/{models[test][0]}/txt_files/2D_v.txt", sep="\s+")
            conv_data["conv_rate"].iloc[0] = np.nan  # Avoid negative diameters
            conv_data["time"] = conv_data["time"]/1.e6

            # New x position for convergence rate circles (adjusted to the left side)
            conv_rate_x_position = -7  # Positioning slightly left of the violin plots

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
                                s=sizes, alpha=0.1, color="grey")  # Grey circles for convergence rate

            xlabel = models[test + "_label"]

            # fig.suptitle(f"{test.capitalize()} Timing")

            ax[0].axhline(y=35, color='k', linestyle='--', linewidth=2)
            ax[1].axhline(y=35, color='k', linestyle='--', linewidth=2)
            ax[0].axvline(x=40, color='grey', linestyle='--', linewidth=1)
            ax[1].axvline(x=40, color='grey', linestyle='--', linewidth=1)
            ax[0].axvline(x=80, color='grey', linestyle='--', linewidth=1)
            ax[1].axvline(x=80, color='grey', linestyle='--', linewidth=1)


            ax[0].set_xlabel(xlabel)
            ax[0].set_ylabel("Time (Myr)")
            ax[1].set_xlabel(xlabel)
            ax[1].set_ylabel("")
            
            plt.tight_layout()  
            plt.savefig(f"{output_dir}/scatter_combined_{test}.png", dpi = 500, format = "png")

           

if __name__ == "__main__":
    main()
