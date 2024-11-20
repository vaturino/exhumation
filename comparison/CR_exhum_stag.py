#!/usr/bin/env python3

# Plot the timing of exhumation, stagnation for velocity models:
################## Convergence rate ####################
##### Time at exhumability ####### Stagnation time ##### 100%
##### Time at exhumability ####### Stagnation time ##### 50%
##### Time at exhumability ####### Stagnation time ##### 0%



import pandas as pd
import numpy as np
import json
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import os
import matplotlib.patches as mpatches


def main():
    colors_tfin = {
        "sed": "cornflowerblue",
        "oc": "#E3B64F",
        "ecl": "#A0C93D",
        "serp": "lightsalmon"
    }


    # Load json file
    with open('models.json') as f:
        models = json.load(f)

    plot_loc = "/home/vturino/PhD/projects/exhumation/plots/single_models"

    tests = ["velocity"] #, "viscosity", "friction", "serpentinization"]
    names = ["velocity_names"]# , "viscosity_names", "friction_names", "serpentinization_names"]


    for idx, test in enumerate(tests):

        pltloc = f"plots/{test}"
        os.makedirs(pltloc, exist_ok=True)

        # # Create the plot
        fig, ax = plt.subplots(4, 2, figsize=(10, 5), gridspec_kw={'height_ratios': [0.25, 1, 1, 1]})

        if test in models:
            model_names = models[names[tests.index(test)]]
            for ind_m, m in enumerate(models[test]):
                text_loc = f"{plot_loc}/{m}/txt_files"
                exhumed_list = pd.read_csv(f"{text_loc}/timing_exhumed_particles.txt", sep="\s+")
                stagnant_list = pd.read_csv(f"{text_loc}/stagnant_particles.txt", sep="\s+")
                cr = pd.read_csv(f"{text_loc}/2D_v.txt", sep="\s+")
                cr["conv_rate"].iloc[0] = np.nan




                labels_font = 7
                ticksize = 6

                edgecolors_exh = "k"
                edgecolors_stag = "k"

            

                alpha = 0.5

        

                # Plot 0: convergence rate
                ax[0, 0].plot(cr["time"]/1e6, cr["conv_rate"], color="grey", label="Convergence rate (cm/yr)", linewidth=1)
                ax[0, 0].axvline(x=35, color="grey", linestyle="--", linewidth = 0.75)
                ax[0, 0].label_outer()  # Only show outer labels and tick labels
                # ax[0, 0].legend(loc='upper right')
                ax[0, 0].set_ylim(0, cr["conv_rate"].max()+1)
                ax[0, 0].set_xlim(0, 50)
                ax[0, 0].tick_params(axis='both', which='major', labelsize=ticksize)

                ax[0, 1].plot(cr["time"]/1e6, cr["conv_rate"], color="grey", label="Convergence rate (cm/yr)", linewidth=1)
                ax[0, 1].axvline(x=35, color="grey", linestyle="--", linewidth = 0.75)
                ax[0, 1].label_outer()  # Only show outer labels and tick labels
                # ax[0, 1].legend(loc='upper right')
                ax[0, 1].set_ylim(0, cr["conv_rate"].max()+1)
                ax[0, 1].set_xlim(0, 50)
                ax[0, 1].tick_params(axis='both', which='major', labelsize=ticksize)


                # Plot 1: Count of particles for each 0.5 Myr for exhumed particles, tfin, do as in first plot
                # ax[ind_m+1, 0].add_patch(plt.Rectangle((0, 0), exhumed_list["tin"].max(), 5000, color='silver', alpha=alpha))
                
                # Create legend handles manually for tfin and tmax using unique lithologies and colors
                unique_lithologies_exhumed = exhumed_list["lithology"].unique()
                handles_tfin_exhumed = [mpatches.Patch(color=colors_tfin[lith], label=lith) for lith in unique_lithologies_exhumed]



                # Histograms for tmax and tfin (Manual layering with plt.bar)
                bin_edges = np.arange(0, 50.5, 1)  # Define bin edges
                width = np.diff(bin_edges)

                sns.histplot(data=exhumed_list, x="tfin", bins=bin_edges, hue="lithology", palette=colors_tfin, ax=ax[ind_m+1, 0], linewidth=0.5, legend = False, element = "step", alpha = alpha)

                # # Plot tfin bars (lower count in front, higher in back)
                # for i in range(len(bin_edges) - 1):
                #     counts = [
                #         (np.histogram(exhumed_list[exhumed_list['lithology'] == lith]['tfin'], bins=bin_edges)[0], colors_tfin[lith]) 
                #         for lith in unique_lithologies_exhumed
                #     ]
                    
                #     # Sort counts by value so that the smaller count bars are drawn in front
                #     sorted_counts = sorted(counts, key=lambda x: x[0][i])
                    
                #     # Plot the smaller count bar (in front)
                #     ax[ind_m+1, 0].bar(bin_edges[i], sorted_counts[0][0][i], width=width[i], color=sorted_counts[0][1],
                #             edgecolor=edgecolors_exh, linewidth=0.5, zorder=3, alpha=1)
                #     # Plot the larger count bar (in back)
                #     ax[ind_m+1, 0].bar(bin_edges[i], sorted_counts[1][0][i], width=width[i], color=sorted_counts[1][1],
                #             edgecolor=edgecolors_exh, linewidth=0.5, zorder=1, alpha=1)

                # Add custom legends for tfin
                # ax[ind_m+1, 0].legend(handles=handles_tfin_exhumed, title="Time at Exhumability", loc=(0.85, 0.85), frameon=True)

                # Set scale, labels, and limits
                ax[ind_m+1, 0].set_yscale("log")
                ax[ind_m+1, 0].label_outer()  # Only show outer labels and tick labels
                # ax[ind_m+1, 0].text(1, 20, " Particles \n Subduction", fontsize=15)
                ax[ind_m+1, 0].set_xlabel("Time (Ma)", fontsize=labels_font)
                ax[ind_m+1, 0].set_ylabel("Particles count (log)", fontsize=labels_font)
                ax[ind_m+1, 0].set_xlim(0, 50)
                ax[ind_m+1, 0].set_ylim(0.9, 5.e3)
                ax[ind_m+1, 0].axvline(x=35, color="grey", linestyle="--", label="35 Ma", linewidth = 0.75)
                ax[ind_m+1, 0].tick_params(axis='both', which='major', labelsize=ticksize)


                # Plot 2: Count of particles for each 0.5 Myr for stagnant particles, tm_kin, tm_dyn, tm_trans: do like for exhumed
                # ax[ind_m+1, 1].add_patch(plt.Rectangle((0, 0), exhumed_list["tin"].max(), 5000, color='silver', alpha=alpha))

                unique_lithologies_stagnant = stagnant_list["lithology"].unique()
                handles_stagnant = [mpatches.Patch(color=colors_tfin[lith], label=lith) for lith in unique_lithologies_stagnant]

                # Define bin edges for tm_kin
                bin_edges = np.arange(0, 50.5, 1)  # Bin edges from 0 to 50 with step 0.5
                width = np.diff(bin_edges)

                sns.histplot(data=stagnant_list, x="tm_kin", bins=bin_edges, hue="lithology", palette=colors_tfin, ax=ax[ind_m+1, 1], alpha=alpha, linewidth=0.5, legend = False, element = "step")
                sns.histplot(data=stagnant_list, x="tm_dyn", bins=bin_edges, hue="lithology", palette=colors_tfin, ax=ax[ind_m+1, 1], alpha=alpha, linewidth=0.5, legend = False, element = "step")
                sns.histplot(data=stagnant_list, x="tm_trans", bins=bin_edges, hue="lithology", palette=colors_tfin, ax=ax[ind_m+1, 1], alpha=alpha, linewidth=0.5, legend = False, element = "step")

                # # Initialize a dictionary to store the histogram counts for each lithology
                # histograms = {}

                # # Calculate histograms for tm_kin for each lithology
                # for lith in unique_lithologies_stagnant:
                #     counts, _ = np.histogram(stagnant_list[stagnant_list['lithology'] == lith]['tm_kin'], bins=bin_edges)
                #     histograms[lith] = counts

                # # Initialize a list to store (count, lithology) pairs
                # sorted_histograms = []

                # # Add all lithology counts to the sorted list
                # for lith in unique_lithologies_stagnant:
                #     sorted_histograms.append((histograms[lith], lith))

                # # Sort the histograms by counts (smaller counts will be first)
                # sorted_histograms.sort(key=lambda x: np.sum(x[0]))  # Sort by the sum of counts for each lithology

                # # Initialize an array for the cumulative counts (to stack the bars)
                # cumulative_counts = np.zeros_like(bin_edges[:-1])

                # # Plot the histograms, ensuring that smaller counts are plotted first
                # for counts, lith in sorted_histograms:
                #     ax[ind_m+1, 1].bar(bin_edges[:-1], counts, width=width, label=lith, color=colors_tfin[lith], 
                #             alpha=1., bottom=cumulative_counts, edgecolor=edgecolors_stag, linewidth = 0.5)  # "bottom" stacks the bars
                #     cumulative_counts += counts  # Update cumulative counts
                # ax[ind_m+1, 1].tick_params(axis='both', which='major', labelsize=ticksize)


                # # Now this but for tm_dyn
                # # Initialize a dictionary to store the histogram counts for each lithology
                # histograms = {}

                # # Calculate histograms for tm_kin for each lithology
                # for lith in unique_lithologies_stagnant:
                #     counts, _ = np.histogram(stagnant_list[stagnant_list['lithology'] == lith]['tm_dyn'], bins=bin_edges)
                #     histograms[lith] = counts

                # # Initialize a list to store (count, lithology) pairs
                # sorted_histograms = []

                # # Add all lithology counts to the sorted list
                # for lith in unique_lithologies_stagnant:
                #     sorted_histograms.append((histograms[lith], lith))

                # # Sort the histograms by counts (smaller counts will be first)
                # sorted_histograms.sort(key=lambda x: np.sum(x[0]))  # Sort by the sum of counts for each lithology

                # # Initialize an array for the cumulative counts (to stack the bars)
                # cumulative_counts = np.zeros_like(bin_edges[:-1])

                # # Plot the histograms, ensuring that smaller counts are plotted first
                # for counts, lith in sorted_histograms:
                #     ax[ind_m+1, 1].bar(bin_edges[:-1], counts, width=width, label=lith, color=colors_tfin[lith], 
                #             alpha=1., bottom=cumulative_counts, edgecolor=edgecolors_stag, linewidth = 0.5)  # "bottom" stacks the bars
                #     cumulative_counts += counts  # Update cumulative counts
                # ax[ind_m+1, 1].tick_params(axis='both', which='major', labelsize=ticksize)




                # # Now this but for tm_trans
                # # Initialize a dictionary to store the histogram counts for each lithology
                # histograms = {}

                # # Calculate histograms for tm_kin for each lithology
                # for lith in unique_lithologies_stagnant:
                #     counts, _ = np.histogram(stagnant_list[stagnant_list['lithology'] == lith]['tm_trans'], bins=bin_edges)
                #     histograms[lith] = counts


                # # Initialize a list to store (count, lithology) pairs
                # sorted_histograms = []

                # # Add all lithology counts to the sorted list
                # for lith in unique_lithologies_stagnant:
                #     sorted_histograms.append((histograms[lith], lith))

                # # Sort the histograms by counts (smaller counts will be first)
                # sorted_histograms.sort(key=lambda x: np.sum(x[0]))  # Sort by the sum of counts for each lithology

                # # Initialize an array for the cumulative counts (to stack the bars)
                # cumulative_counts = np.zeros_like(bin_edges[:-1])

                # # Plot the histograms, ensuring that smaller counts are plotted first
                # for counts, lith in sorted_histograms:
                #     ax[ind_m+1, 1].bar(bin_edges[:-1], counts, width=width, label=lith, color=colors_tfin[lith], 
                #             alpha=1., bottom=cumulative_counts, edgecolor=edgecolors_stag, linewidth = 0.5)  # "bottom" stacks the bars
                #     cumulative_counts += counts  # Update cumulative counts
                # ax[ind_m+1, 1].tick_params(axis='both', which='major', labelsize=ticksize)





                # Add a title, labels, and legend
                # ax[ind_m+1, 1].text(1, 20, " Particles \n Subduction", fontsize=15)
                # ax[ind_m+1, 1].legend(handles=handles_stagnant, title="Lithologies")
                ax[ind_m+1, 1].set_xlim(0, 50)
                ax[ind_m+1, 1].set_ylim(0.9, 5.e3)
                ax[ind_m+1, 1].set_yscale("log")
                ax[ind_m+1, 1].axvline(x=35, color="grey", linestyle="--", label="35 Ma", linewidth = 0.75)
                ax[ind_m+1, 1].tick_params(axis='both', which='major', labelsize=ticksize)



            

            # Add no vertical space between the plots
            plt.subplots_adjust(hspace=0)

            plt.savefig(f"{pltloc}/CR_exhumed_stagnant.pdf", dpi=500)
            plt.close()





if __name__ == "__main__":
    main()

