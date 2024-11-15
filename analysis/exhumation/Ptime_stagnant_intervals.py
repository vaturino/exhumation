#!/usr/bin/env python3

import pandas as pd
import numpy as np
import json
import matplotlib.pyplot as plt
from pathlib import Path
import argparse
import os
import logging

############### FUNCTIONS ####################

def setup_logging():
    """Setup logging configuration."""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

def load_json_config(json_file, json_loc):
    """Load configuration from JSON file."""
    json_path = Path(json_loc) / json_file
    if not json_path.exists():
        logging.error(f"JSON file not found: {json_path}")
        raise FileNotFoundError(f"JSON file not found: {json_path}")
    
    with open(json_path, 'r') as file:
        return json.load(file)

def create_directory(path):
    """Create directory if it doesn't exist."""
    os.makedirs(path, exist_ok=True)

def validate_columns(df, required_columns):
    """Check if required columns are in the DataFrame."""
    missing = [col for col in required_columns if col not in df.columns]
    if missing:
        logging.warning(f"Missing columns in the data: {missing}")
        return False
    return True

############### MAIN ####################

def main():
    setup_logging()

    # Define colors for lithologies
    colors_tmax = {
        "sed": "mediumblue",  # Color for sediment
        "oc": "#B06D1A",  # Color for organic
        "ecl": "#45701C",  # Color for eclogite
        "serp": "brown"  # Color for serpentinite
    }

    colors_tfin = {
        "sed": "cornflowerblue",  # Lighter color for error shading
        "oc": "#E3B64F",
        "ecl": "#A0C93D",
        "serp": "lightsalmon"
    }

    # Alpha values for controlling transparency
    alpha_dist = 0.1  # Alpha for main polygon (distribution)
    alpha_err = 0.15  # Alpha for error shading

    parser = argparse.ArgumentParser(description='Script to process and analyze particle pressure bins by lithology.')
    parser.add_argument('json_file', help='JSON file with model name, time at the end of computation, folder where to save the plots, legend')
    args = parser.parse_args()

    json_loc = '/home/vturino/PhD/projects/exhumation/pyInput/'
    configs = load_json_config(args.json_file, json_loc)
    model_name = configs["models"][0]

    plot_loc = f"/home/vturino/PhD/projects/exhumation/plots/single_models/{model_name}"
    txt_loc = f"{plot_loc}/txt_files"
    create_directory(txt_loc)

    data_file = f"{txt_loc}/stagnant_particles.txt"
    if not Path(data_file).exists():
        logging.error(f"Data file not found: {data_file}")
        return

    # Read the data from the CSV file
    data = pd.read_csv(data_file, sep=r"\s+")
    cr = pd.read_csv(f"{txt_loc}/2D_v.txt", sep=r"\s+")
    cr.iloc[0] = np.nan
    cr = cr.dropna()

    if data.empty:
        logging.warning(f"Data file is empty: {data_file}")
        return

    # Round pressures to 0.1 increments
    for i in ["kin", "dyn", "trans"]:
        data[f"Pm_{i}"] = np.round(data[f"Pm_{i}"], 1)
        data[f"ti_{i}"] = np.round(data[f"ti_{i}"], 1)
        data[f"tf_{i}"] = np.round(data[f"tf_{i}"], 1)

    # Group by lithology and pressure and calculate mean, min, and max for ti_i and tf_i
    lithologies = data["lithology"].unique()  # Get unique lithologies

    # Plot
    # f1, a1 = plt.subplots(3, 1, figsize=(17, 15), gridspec_kw={'height_ratios': [0.5, 1.5, 1.5]})
    f1,a1 = plt.subplots(2, 1, figsize=(12, 8), gridspec_kw={'height_ratios': [0.25, 1]})

    # Plot convergence rate in the first subplot
    a1[0].plot(cr['time'] / 1e6, cr['conv_rate'], color='grey', label='Convergence rate', linewidth=2)
    a1[0].label_outer()
    a1[0].set_xlim(0, 50)
    a1[0].set_ylim(0, max(cr['conv_rate']) + 0.2)
    a1[0].axvline(x=35., color='grey', linewidth=1, linestyle="--")
    a1[0].legend()
    a1[0].label_outer()

    # Define a flag to track whether the legend has been added for each lithology
    for lithology in lithologies:
        # Filter data by lithology
        lithology_data = data[data["lithology"] == lithology]

        for i in ["kin", "dyn", "trans"]:
            # Group by pressure and calculate mean, min, and max for ti_i and tf_i
            grouped = (
                lithology_data.groupby(f"Pm_{i}")
                    .agg(
                        ti_mean=(f"ti_{i}", "mean"),
                        tf_mean=(f"tf_{i}", "mean"),
                        ti_min=(f"ti_{i}", "min"),
                        ti_max=(f"ti_{i}", "max"),
                        tf_min=(f"tf_{i}", "min"),
                        tf_max=(f"tf_{i}", "max"),
                        tm_mean=(f"tm_{i}", "mean"),
                        sediment_count=(f"Pm_{i}", "count")
                    )
                    .reset_index()
            )

            # Check if the grouped DataFrame is empty
            if grouped.empty:
                logging.warning(f"Grouped data for {lithology} {i} is empty.")
                continue

            # Calculate lower and upper bounds based on min and max (variation range)
            grouped["ti_lower"] = grouped["ti_min"]
            grouped["ti_upper"] = grouped["ti_max"]
            grouped["tf_lower"] = grouped["tf_min"]
            grouped["tf_upper"] = grouped["tf_max"]

            # Extracting x and y coordinates for the polygon
            x_bottom = grouped["ti_mean"]  # Using the mean for the bottom contour (ti)
            x_top = grouped["tf_mean"]  # Using the mean for the top contour (tf)
            x_bottom_lower = grouped["ti_lower"]
            x_bottom_upper = grouped["ti_upper"]
            x_top_lower = grouped["tf_lower"]
            x_top_upper = grouped["tf_upper"]
            y = grouped[f"Pm_{i}"]

            # Find the index of the maximum pressure
            max_pressure_idx = grouped[f"Pm_{i}"].idxmax()

            # Get the ti_mean and tf_mean values for the highest pressure
            ti_max_pressure = grouped["ti_mean"].iloc[max_pressure_idx]
            tf_max_pressure = grouped["tf_mean"].iloc[max_pressure_idx]
            max_pressure = grouped[f"Pm_{i}"].iloc[max_pressure_idx]

            # Check if there is data available to plot the line
            if len(grouped) > 0:
                # Plotting the custom bottom line you requested (connecting the mean ti and tf values)
                a1[1].plot([grouped["ti_mean"].iloc[0], grouped["tf_mean"].iloc[0]], 
                        [grouped[f"Pm_{i}"].iloc[0], grouped[f"Pm_{i}"].iloc[0]], 
                        color=colors_tmax[lithology], 
                        linewidth=1, linestyle='--')

                # Plotting the custom top line (connecting the mean tf and ti values at the maximum pressure)
                a1[1].plot([ti_max_pressure, tf_max_pressure], 
                        [max_pressure, max_pressure], 
                        color=colors_tmax[lithology], 
                        linewidth=1, linestyle='--')

                # Fill polygon for the main region (using the lithology color)
                a1[1].fill_betweenx(y, x_bottom, x_top, color=colors_tmax[lithology], alpha=alpha_dist)

                # Add shading for variation range errors with the lighter lithology color
                a1[1].fill_betweenx(y, x_bottom_lower, x_bottom_upper, color=colors_tfin[lithology], alpha=alpha_err)
                a1[1].fill_betweenx(y, x_top_lower, x_top_upper, color=colors_tfin[lithology], alpha=alpha_err)

                # Plot lines for all polygon sides using the lithology color
                a1[1].plot(x_bottom, y, color=colors_tmax[lithology], linewidth=1, linestyle= "--")
                a1[1].plot(x_top, y, color=colors_tmax[lithology], linewidth=1, linestyle= "--")

                # Plot scatter with size proportional to sediment_count
                a1[1].scatter(grouped["tm_mean"], grouped[f"Pm_{i}"], color=colors_tmax[lithology], 
                            s=grouped["sediment_count"], alpha=0.7, edgecolors='black')

    # Labels and title for the main plot
    a1[1].set_xlabel("Time (Myr)", fontsize=15)
    a1[1].set_ylabel("Pressure (GPa)", fontsize=15)
    a1[1].set_xlim(0, 50)
    a1[1].set_ylim(0, 2.6)
    a1[1].axvline(x=35., color='grey', linewidth=1, linestyle="--")
    # Add legend for the scatter plot
    handles, labels = a1[1].get_legend_handles_labels()
    scatter_legend = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=colors_tmax[lith], markersize=10, label=lith) for lith in lithologies]
    a1[1].legend(handles=scatter_legend, title="Lithology", loc='upper left')

    f1.subplots_adjust(hspace=0.)
    plt.savefig(f"{plot_loc}/stagnation_distribution.png", dpi=500)
    plt.close()

    # if 'tm_mean' not in data.columns:
    #     # Example: If 'tm_mean' is the mean of 'ti_mean' and 'tf_mean'
    #     data['tm_mean'] = (data['ti_mean'] + data['tf_mean']) / 2  # Modify this logic as needed

    # # Now, proceed with grouping by tm_mean and calculating the particle count
    # grouped_by_tm = data.groupby("tm_mean")["sediment_count"].sum().reset_index(name='count')

    # # Plot histogram for particle count against tm_mean
    # for lithology in lithologies:
    #     # Filter by lithology
    #     lithology_data = grouped_by_tm[grouped_by_tm["lithology"] == lithology]

    #     # Plot the histogram for this lithology in a1[2]
    #     a1[2].bar(lithology_data["tm_mean"], lithology_data["count"], width=0.1, color=colors_tmax[lithology], label=lithology, alpha=0.7)

    # # Set labels and title for the histogram
    # a1[2].set_xlabel("Time (tm_mean)")
    # a1[2].set_ylabel("Particle Count (sediment_count)")
    # a1[2].set_title("Histogram of Particle Count by Time (tm_mean) and Lithology")

    # # Add legend
    # a1[2].legend(title="Lithology")





if __name__ == "__main__":
    main()