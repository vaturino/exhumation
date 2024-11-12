#! /usr/bin/python3

from matplotlib.cm import get_cmap
import pandas as pd
import numpy as np
import json
from tqdm import tqdm
import matplotlib.cm as cm
import matplotlib as mpl
import sys, os, subprocess
from matplotlib.pyplot import copper, show, xlabel, ylabel
import matplotlib.pyplot as plt
import argparse
from matplotlib.gridspec import GridSpec
from matplotlib import colors
import math as math
from scipy.signal import savgol_filter 
import seaborn as sns
from matplotlib.axes._axes import _log as matplotlib_axes_logger
matplotlib_axes_logger.setLevel('ERROR')
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
pd.options.mode.chained_assignment = None  # default='warn'
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

path = str(Path(Path(__file__).parent.absolute()).parent.absolute())
sys.path.insert(0, path)
from libraries.functions import *
from libraries.particles import *
from libraries.exhumation import *  
import plotly.express as px


def load_data(id, txt_loc):
    data = pd.read_csv(f"{txt_loc}/PT/pt_part_{id}.txt", sep="\s+")
    data["time"] = data["time"] / 2
    data["ts"] = data["time"] * 2
    return data



def main():
    parser = argparse.ArgumentParser(description='Script that gets some models and gives the kinematic indicators')
    parser.add_argument('json_file', help='json file with model name, time at the end of computation, folder where to save the plots, legend')
    args = parser.parse_args()

    json_loc = '/home/vturino/PhD/projects/exhumation/pyInput/'

    # Read the json file
    with open(f"{json_loc}{args.json_file}") as json_file:
        configs = json.load(json_file)
    m = configs["models"][0]

    # Create the folders to save the plots
    plot_loc = f"/home/vturino/PhD/projects/exhumation/plots/single_models/{m}"
    txt_loc = f'{plot_loc}/txt_files'
    if not os.path.exists(txt_loc):
        os.mkdir(txt_loc)

    sloc = f"{plot_loc}/stagnant_times"
    if not os.path.exists(sloc):
        os.mkdir(sloc)

    colors_tmax = {
        "sed": "mediumblue",
        "oc": "#B06D1A",
        "ecl": "#45701C",
        "serp": "brown"
    }


    
    part = pd.read_csv(f"{txt_loc}/stagnant_particles.txt", sep="\s+")
    particles = pd.DataFrame(columns=["id", "lithology", "Pm", "tm", "time_interval", "ti", "tf"], index=range(len(part)))

    for i in range(len(part)):
    # for i in [18]:
        id = part["id"].iloc[i]
        # print(id)
        data = load_data(id, txt_loc)

        grad_thresh = 0.01
        time_thresh = 0.5
        data = data[data["time"] > 10]
        data['gradient'] = np.gradient(data["Plith"], data["time"])

        # find time intervals where the gradient is below the threshold
        lowgrad = data[abs(data["gradient"]) < grad_thresh].reset_index(drop=True)
        
        # Calculate the duration of each interval (time difference between consecutive rows)
        lowgrad["duration"] = lowgrad["time"].diff()

        # Initialize the time_bin column
        lowgrad['time_bin'] = None

        start_time = None
        for i, row in lowgrad.iterrows():
            # If duration == time_thresh, we're still in the same interval
            if row['duration'] == time_thresh:
                if start_time is None:
                    start_time = row['time']  # Mark the start time of the interval
                end_time = row['time']  # Update the end time (current row's time)

                # If we reach the last row of consecutive 0.5 intervals or the last row in the bin, assign the time_bin
                if i == len(lowgrad) - 1 or lowgrad.at[i+1, 'duration'] != time_thresh:
                    # Assign the time bin to all rows in the current interval
                    for j in range(i, -1, -1):
                        if lowgrad.at[j, 'duration'] == time_thresh:
                            if start_time != end_time:  
                                lowgrad.at[j, 'time_bin'] = f"[{start_time}, {end_time})"
                        else:
                            break
                    start_time = None  # Reset start_time after assigning the time bin
                else:
                    # Continue updating end_time when in a consecutive 0.5 interval
                    continue
            else:
                start_time = None  # Reset when a non-0.5 duration is encountered

        # Group by 'time_bin' and compute the average of 'Plith' and 'time' for each interval
        avg_values = lowgrad.groupby('time_bin').agg({'Plith': 'mean', 'time': 'mean'}).reset_index()
        lowgrad["Pm"] = np.nan
        lowgrad["tm"] = np.nan

        for tint in avg_values["time_bin"]:
            # Get the average 'Plith' value for the current time bin
            lowgrad["Pm"] = np.where(lowgrad["time_bin"] == tint, avg_values[avg_values["time_bin"] == tint]["Plith"].values[0], lowgrad["Pm"])
            lowgrad["tm"] = np.where(lowgrad["time_bin"] == tint, avg_values[avg_values["time_bin"] == tint]["time"].values[0], lowgrad["tm"])
     
        
        if i % 1 == 0:
            fig = plt.figure(figsize=(8, 6))
            plt.plot(lowgrad["time"], lowgrad["Pm"], color="green", linewidth = 2, zorder = 5)
            plt.plot(data["time"], data["Plith"], color="grey")
            plt.scatter(lowgrad["tm"], lowgrad["Pm"], color="red", zorder = 10)
            plt.xlabel("Time (Myr)")
            plt.ylabel("Pressure (GPa)")
            plt.title(f"Particle {id} - Lithology: {part['lithology'].iloc[i]}")
            plt.ylim(0,2.5)
            plt.savefig(f"{sloc}/stagnant_times_{id}.png")
            plt.close()

    # particles.dropna(subset=["id"], inplace=True)
    # print("Number of stagnant particles: ", len(particles))

    # float_cols = particles.columns.difference(["id", "lithology"])
    # particles[float_cols] = particles[float_cols].astype(float)
    # particles.to_csv(f"{txt_loc}/stagnant_times.txt", sep="\t", index=False, header=True, float_format='%.2f')





if __name__ == "__main__":
    main()


#    # Find the lowest, mid and highest values of the time interval for each lithology,
#    # at the end I want a dataframe like this:
#     # lithology | time_interval | avg_tm | avg_tin | avg_tf | count
#     # --------------------------------------------------------------
#     # sed       | 10            | 20     | 5       | 15     | 100
#     # sed       | 15            | 20     | 5       | 25     | 50
#     # sed       | 20            | 20     | 5       | 30     | 25
#     # ecl       | 25            | 20     | 5       | 35     | 10
#     # ecl       | 30            | 20     | 5       | 40     | 5
#     # ecl       | 35            | 20     | 5       | 45     | 5
#     # oc        | 5             | 20     | 5       | 10     | 100
#     # oc        | 10            | 20     | 5       | 15     | 100
#     # oc        | 30            | 20     | 5       | 40     | 5

#     # Define the number of bins
#     num_bins = 5

#     # Initialize an empty DataFrame to store the results
#     pressure_bins = pd.DataFrame()

#     # Loop through each unique lithology
#     for lithology in particles["lithology"].unique():
#         # Filter particles by lithology
#         lithology_particles = particles[particles["lithology"] == lithology]
        
#         # Find equal intervals of Pm for the current lithology
#         pressures = lithology_particles["Pm"]
#         bins = pd.cut(pressures, num_bins).unique()
        
#         # Create a DataFrame for the current lithology and append to the result
#         lithology_bins = pd.DataFrame({
#             "lithology": [lithology] * len(bins),
#             "pressure_bin": bins
#         })
#         pressure_bins = pd.concat([pressure_bins, lithology_bins], ignore_index=True)

#     # add column with count of values for each bin
#     pressure_bins["count"] = pressure_bins.apply(
#         lambda row: len(particles[(particles["lithology"] == row["lithology"]) & (pd.cut(particles["Pm"], bins=[row["pressure_bin"]]).notna())]),
#         axis=1
#     )

#     # add columns with average value of time_interval, ti and tf for each bin
#     for col in ["time_interval", "ti", "tf"]:
#         pressure_bins[f"avg_{col}"] = pressure_bins.apply(
#             lambda row: particles[(particles["lithology"] == row["lithology"]) & (pd.cut(particles["Pm"], bins=[row["pressure_bin"]]).notna())][col].mean(),
#             axis=1
#         )

#     # tm is the average of time_interval + ti
#     pressure_bins["avg_tm"] = pressure_bins["avg_ti"] + pressure_bins["avg_time_interval"]

#     # Plot the results as scatter + error bar: x-axis = pressure, y-axis = avg_tm, color = lithology, error = avg_time_interval
#     fig = plt.figure(figsize=(10, 6))
    




#     print(pressure_bins)
#     exit()
    

    # # Group by lithology and time_interval, then compute the average `tm` and `tin`
    # grouped = particles.groupby(['lithology', 'time_interval']).agg(
    #     avg_tm=('tm', 'mean'),
    #     avg_tin=('ti', 'mean'),
    #     avg_tf=('tf', 'mean'),
    #     count=('time_interval', 'size')
    # ).reset_index()

    # # Get the top 5 most frequent intervals for each lithology
    # top_intervals = (
    #     grouped.sort_values(['lithology', 'count', 'time_interval'], ascending=[True, False, False])
    #     .groupby('lithology')
    #     .head(5)
    # )

    # # Plotting
    # plt.figure(figsize=(10, 6))

    # unique_lithologies = top_intervals['lithology'].unique()
    # # Create a list to hold custom legend handles for size
    # legend_handles = []

    # # Loop through each lithology and apply unique offsets
    # for lith in unique_lithologies:
    #     lith_group = top_intervals[top_intervals['lithology'] == lith]
    #     n = len(lith_group)
        
    #     # Generate offsets for each entry in the group
    #     offsets = np.linspace(-0.2, 0.2, n)
        
    #     for offset, (_, row) in zip(offsets, lith_group.iterrows()):
    #         # Apply offset to avoid overlap
    #         lithology_position = unique_lithologies.tolist().index(row['lithology']) + offset
    #         tmed = row["time_interval"] / 2 + row["avg_tin"]
    #         color = colors_tmax.get(row['lithology'], "gray")
            
    #         # Plot the time interval as a line from avg_tin to avg_tin + time_interval
    #         plt.plot(
    #             [lithology_position, lithology_position],
    #             [row['avg_tin'] , row['avg_tin'] + row["time_interval"]],
    #             color=color, alpha=1.
    #         )
            
    #         # Mark avg_tm as a point within the time interval
    #         plt.scatter(lithology_position, tmed, 
    #                     color=color, s = row["count"]*2, 
    #                     zorder=5, 
    #                     label='avg_tm'
    #                     ) # if (lith == unique_lithologies[0] and offset == offsets[0]) else "")
            
            

    # # Customizing the plot
    # plt.xticks(ticks=np.arange(len(unique_lithologies)), labels=unique_lithologies)
    # plt.xlabel("Lithology")
    # plt.ylabel("Time")
    # plt.title("Average time interval per lithology")
    # # Add custom legend for size of the points
    # # plt.legend(title="Count of points", loc='upper right')
    # # legend_handles = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=colors_tmax[lith], markersize=8, label=lith)
    # #                 for lith in unique_lithologies]
    # # plt.legend(handles=legend_handles, title="Lithology")

    # plt.ylim(0, 50)

    # plt.savefig(f"{plot_loc}/stagnant_times.png")
    # plt.close()

    # print(particles[(particles["lithology"] == "sed") & (particles["tf"] < 49.5)])





# def compute_time_intervals(data, Pthresh, min_delta_time=10):
#     intervals = []
#     start, end = None, None
    
#     for i in range(len(data) - 1):
#         if data["time"].iloc[i] > 10 and abs(data["Plith"].iloc[i] - data["Plith"].iloc[i + 1]) <= Pthresh:
#             if start is None:
#                 start = i
#             end = i + 1
#         else:
#             if start is not None and (data["time"].iloc[end] - data["time"].iloc[start]) >= min_delta_time:
#                 Pm = data["Plith"].iloc[start:end + 1].mean()
#                 interval = {
#                     "start_time": data["time"].iloc[start],
#                     "end_time": data["time"].iloc[end],
#                     "Pm": Pm,
#                     "duration": data["time"].iloc[end] - data["time"].iloc[start]
#                 }
#                 intervals.append(interval)
#             start, end = None, None

#     if start is not None and (data["time"].iloc[end] - data["time"].iloc[start]) >= min_delta_time:
#         Pm = data["Plith"].iloc[start:end + 1].mean()
#         interval = {
#             "start_time": data["time"].iloc[start],
#             "end_time": data["time"].iloc[end],
#             "Pm": Pm,
#             "duration": data["time"].iloc[end] - data["time"].iloc[start]
#         }
#         intervals.append(interval)

#     return intervals


        # intervals = compute_time_intervals(data, thresh)
        # rows_to_add = []

        # for interval in intervals:
        #     if interval["duration"] >= 10:
        #         rows_to_add.append({
        #             "id": id,
        #             "lithology": part["lithology"].iloc[i],
        #             "Pm": interval["Pm"],
        #             "tm": (interval["start_time"] + interval["end_time"]) / 2,
        #             "time_interval": interval["duration"]
        #         })

        # if rows_to_add:
        #     particles = pd.concat([particles, pd.DataFrame(rows_to_add)], ignore_index=True)