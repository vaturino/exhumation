#! /usr/bin/python3

import pandas as pd
import seaborn as sns
import json
import joypy



import matplotlib.pyplot as plt

def main():
    # Load json file
    with open('models.json') as f:
        models = json.load(f)

    plot_loc = "/home/vturino/PhD/projects/exhumation/plots/single_models"
    tests = ["velocity", "viscosity", "friction", "serpentinization"]
    names = ["velocity_names", "viscosity_names", "friction_names", "serpentinization_names"]

    # Set colorblind-friendly palette
    palette = sns.color_palette("colorblind")
    
    # Define different line styles for lithology values
    linestyles = ['-', '--', '-.', ':']  # Four different styles
    
    for test in tests:
        if test in models:
            model_names = models[names[tests.index(test)]]
            
            # Create figure and subplots once for all models
            fig, ax = plt.subplots(2, 3, figsize=(15, 10))
            fig.suptitle(f"Timing of exhumation and stagnation as a function of {test}")

            # Loop through each model to plot their respective data
            for ind_m, m in enumerate(models[test]):
                text_loc = f"{plot_loc}/{m}/txt_files"
                
                try:
                    exhumed = pd.read_csv(f"{text_loc}/timing_exhumed_particles.txt", sep="\s+")
                    stagnant = pd.read_csv(f"{text_loc}/timing_stagnant_particles.txt", sep="\s+")
                    cr = pd.read_csv(f"{text_loc}/2D_v.txt", sep="\s+")
                except FileNotFoundError:
                    print(f"Files not found for model {m}. Skipping.")
                    continue

                smooth = 1.5
                # normalize time by convergence rate time
                exhumed["tin"] = exhumed["tin"]/cr["time"].iloc[-1]*1e6
                exhumed["tmax"] = exhumed["tmax"]/cr["time"].iloc[-1]*1e6
                exhumed["tfin"] = exhumed["tfin"]/cr["time"].iloc[-1]*1e6
                
                stagnant["tin"] = stagnant["tin"]/cr["time"].iloc[-1]*1e6
                stagnant["tmax"] = stagnant["tmax"]/cr["time"].iloc[-1]*1e6
                stagnant["tfin"] = stagnant["tfin"]/cr["time"].iloc[-1]*1e6


                # print(exhumed[['tin', 'tmax', 'tfin']].describe())
                # print(stagnant[['tin', 'tmax', 'tfin']].describe())
                # exit()

               
                ### Exhumed particles ###
                # Get unique lithology values
                unique_lithologies = exhumed["lithology"].unique()

                # Get the color for this model based on its index (ind_m)
                color = palette[ind_m % len(palette)]



                # Plot kdeplot for each lithology for tin
                for i, lith in enumerate(unique_lithologies):
                    linestyle = linestyles[i % len(linestyles)]  # Cycle through linestyles
                    # ax[0,0].hist(exhumed[exhumed["lithology"] == lith]["tin"], bins=30, density=True, alpha=0.5)

                    sns.kdeplot(
                        exhumed[exhumed["lithology"] == lith]["tin"], 
                        ax=ax[0,0], 
                        label=f'{model_names[ind_m]} - {lith}', 
                        linestyle=linestyle,
                        color=color,  
                        common_norm=True,  # Normalize the distribution
                        bw_adjust=smooth
                    )

                # Plot kdeplot for each lithology for tmax
                for i, lith in enumerate(unique_lithologies):
                    linestyle = linestyles[i % len(linestyles)]  # Cycle through linestyles
                    sns.kdeplot(
                        exhumed[exhumed["lithology"] == lith]["tmax"], 
                        ax=ax[0,1], 
                        label=f'{model_names[ind_m]} - {lith}', 
                        linestyle=linestyle,
                        color=color,  
                        common_norm=True,  # Normalize the distribution
                        bw_adjust=smooth
                    )

                # Plot kdeplot for each lithology for tfin
                for i, lith in enumerate(unique_lithologies):
                    linestyle = linestyles[i % len(linestyles)]
                    sns.kdeplot(
                        exhumed[exhumed["lithology"] == lith]["tfin"], 
                        ax=ax[0,2], 
                        label=f'{model_names[ind_m]} - {lith}', 
                        linestyle=linestyle,
                        color=color, 
                        common_norm=True,  # Normalize the distribution
                        bw_adjust=smooth
                    )

                ### Stagnant particles ###
                # Get unique lithology values
                unique_lithologies = stagnant["lithology"].unique()

                # Plot kdeplot for each lithology for tin
                for i, lith in enumerate(unique_lithologies):
                    linestyle = linestyles[i % len(linestyles)]
                    sns.kdeplot(
                        stagnant[stagnant["lithology"] == lith]["tin"], 
                        ax=ax[1,0], 
                        label=f'{model_names[ind_m]} - {lith}', 
                        linestyle=linestyle,
                        color=color, 
                        common_norm=True,  # Normalize the distribution
                        bw_adjust=smooth
                    )
                
                # Plot kdeplot for each lithology for tmax
                for i, lith in enumerate(unique_lithologies):
                    linestyle = linestyles[i % len(linestyles)]
                    sns.kdeplot(
                        stagnant[stagnant["lithology"] == lith]["tmax"], 
                        ax=ax[1,1], 
                        label=f'{model_names[ind_m]} - {lith}', 
                        linestyle=linestyle,
                        color=color,  
                        common_norm=True,  # Normalize the distribution
                        bw_adjust=smooth
                    )

                # Plot kdeplot for each lithology for tfin
                for i, lith in enumerate(unique_lithologies):
                    linestyle = linestyles[i % len(linestyles)]
                    sns.kdeplot(
                        stagnant[stagnant["lithology"] == lith]["tfin"], 
                        ax=ax[1,2], 
                        label=f'{model_names[ind_m]} - {lith}', 
                        linestyle=linestyle,
                        color=color, 
                        common_norm=True,  # Normalize the distribution
                        bw_adjust=smooth
                    )
            
            # Finalize ax[0,0] by adding labels and legend
            ax[0,0].set_title("Exhumed particles")
            ax[0,0].set_xlabel("Subduction time (Ma)")
            ax[0,0].set_ylabel("Density")
            ax[0,0].legend(title="Model - Lithology")
            

            ax[0,1].set_xlabel("Time at peak pressure (Ma)")
            ax[0,1].set_ylabel("Density")

            ax[0,2].set_xlabel("Time at exhumation (Ma)")
            ax[0,2].set_ylabel("Density")

            # Finalize axes for stagnant particles
            ax[1,0].set_title("Stagnant particles")
            ax[1,0].set_xlabel("Subduction time (Ma)")
            ax[1,0].set_ylabel("Density")
            

            ax[1,1].set_xlabel("Time at peak pressure (Ma)")
            ax[1,1].set_ylabel("Density")

            ax[1,2].set_xlabel("Time at stagnation (Ma)")
            ax[1,2].set_ylabel("Density")  
            ax[1,2].legend(title="Model - Lithology", loc='upper left') 

            # Normalize y-axis to range from 0 to 1
            for i in range(2):
                for j in range(3):
                    ax[i,j].set_xlim(0, 1)

            # Save or display the plot
            plt.tight_layout()
            plt.subplots_adjust(top=0.9)  # Adjust layout to fit title
            output_file = f"combined/timing/timing_density_plots_{test}.png"
            plt.savefig(output_file)
            plt.close()

if __name__ == '__main__':
    main()
