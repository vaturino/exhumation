import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tqdm import tqdm
import json
import argparse
import os

def main():
    parser = argparse.ArgumentParser(description='Script that gets some models and gives the kinematic indicators')
    parser.add_argument('json_file', help='json file with model name, time at the end of computation, folder where to save the plots, legend')
    args = parser.parse_args()

    json_loc = '/home/vturino/PhD/projects/exhumation/pyInput/'

    with open(f"{json_loc}{args.json_file}") as json_file:
        configs = json.load(json_file)
    m = configs["models"][0]

    plot_loc = f"/home/vturino/PhD/projects/exhumation/plots/single_models/{m}"
    txt_loc = f'{plot_loc}/txt_files'
    if not os.path.exists(txt_loc):
        os.mkdir(txt_loc)

    eloc = f"{plot_loc}/exhumed_Ptime"
    if not os.path.exists(eloc):
        os.mkdir(eloc)

    # Load the data
    data = pd.read_csv(f"{txt_loc}/exhumed_particles.txt", sep="\s+")

    final_data = pd.DataFrame(columns=['id', 'tin', 'tfin', 'tmid', 'Pin', 'Pfin','Pmid', 'time_interval', 'lithology'])

    for p in tqdm(range(len(data))):
        id = data["id"].iloc[p]
        particle = pd.read_csv(f"{txt_loc}/PT/pt_part_{id}.txt", sep="\s+")

        # Parameters
        Pmax = particle['Plith'].max()
        Pthresh = round(0.75 * Pmax, 3)
        grad_thresh = 0.05
        tmax = particle['time'][particle['Plith'].idxmax()]
        post_tmax_indices = particle['time'] > tmax

        # Get the pressures and times after tmax
        post_tmax_times = particle['time'][post_tmax_indices]
        post_tmax_Plith = particle['Plith'][post_tmax_indices]

        # Find the time when Plith first drops below 25% of max Plith
        below_thresh_indices = post_tmax_Plith <= Pthresh

        if below_thresh_indices.any():
            # Get the first time index where Plith is below 0.75 * Pmax
            tinit = post_tmax_times[below_thresh_indices.idxmax()]
        else:
            tinit = None

        # Process intervals where Plith >= 0.75 * Pmax after tmax
        if tinit is not None:
            # Find when Plith is >= 0.75 * Pmax after tmax
            above_thresh_indices = post_tmax_Plith >= Pthresh

            if above_thresh_indices.any():
                # Get the first time where Plith exceeds or equals 0.75 * Pmax
                tstart = post_tmax_times[above_thresh_indices.idxmax()]
                post_tmax_times_valid = post_tmax_times[post_tmax_times >= tstart]
                post_tmax_Plith_valid = post_tmax_Plith[post_tmax_times >= tstart]

                # Compute the gradient of pressure with respect to time
                pressure_gradient = np.gradient(post_tmax_Plith_valid, post_tmax_times_valid)

                # Filter the gradient and times based on the condition (gradient < 1e-2)
                valid_indices = (np.abs(pressure_gradient) < grad_thresh)
                stabilizing_index = np.where(valid_indices)[0]

                if stabilizing_index.size > 0:
                    t_stabilized = post_tmax_times_valid.iloc[stabilizing_index[0]]
                else:
                    t_stabilized = None

                # Calculate the time when the gradient exceeds the threshold again
                valid_indices_above_threshold = (np.abs(pressure_gradient) > grad_thresh)
                first_above_threshold_index = np.where(valid_indices_above_threshold)[0]

                if first_above_threshold_index.size > 0:
                    t_threshold_breach = post_tmax_times_valid.iloc[first_above_threshold_index[0]]
                else:
                    t_threshold_breach = post_tmax_times_valid.max()

                # Calculate middle time
                if t_stabilized and t_threshold_breach:
                    t_middle = (t_stabilized + t_threshold_breach) / 2
                    p_middle = particle['Plith'][particle['time'] == t_middle].values[0]
                    p_init = particle['Plith'][particle['time'] == t_stabilized].values[0]
                    p_breach = particle['Plith'][particle['time'] == t_threshold_breach].values[0]
                    time_interval = t_threshold_breach - t_stabilized

                    # Save the results in final_data DataFrame
                    final_data.loc[p] = [id, tinit/2, t_threshold_breach/2, t_middle/2, p_init, p_breach, p_middle, time_interval/2, data['lithology'].iloc[p]]
                
                # Plotting the pressure over time for validation
                plt.plot(particle["time"]/2, particle["Plith"], label='Pressure', color='blue')

                # Adding labels and title
                if p % 10 == 0:
                    plt.scatter(t_stabilized/2, p_init, color='red', label='Stabilized Time')
                    plt.scatter(t_threshold_breach/2, p_breach, color='orange', label='Threshold Breach Time')
                    plt.scatter(t_middle/2, p_middle, color='green', label='Middle Time')
                    plt.xlabel('Time (Myr)')
                    plt.ylabel('Plith')
                    plt.title(f'Pressure vs Time for Particle {id} - Lithology: {data["lithology"].iloc[p]}')
                    plt.legend()
                    plt.tight_layout()
                    plt.savefig(f"{eloc}/pressure_vs_time_particle_{id}.png")
                plt.clf()

            else:
                print(f"No valid interval for particle {id} where Plith >= 0.75 * Pmax")

        else:
            print(f"Particle {id} - tinit is None, no data to plot.")
            continue

    # Save the final data to a file
    final_data.to_csv(f"{txt_loc}/exhumed_particles_timing.txt", sep="\t", index=False)

if __name__ == "__main__":
    main()
