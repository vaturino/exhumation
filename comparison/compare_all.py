#! /usr/bin/python3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def main():
    # Load the data from a spreadsheet
    comparison = pd.read_excel('comparison.xlsx')

    # Subset the data by test types
    viscosity = comparison[(comparison["Test"] == "Viscosity") | (comparison["Test"] == "Reference")]
    friction = comparison[(comparison["Test"] == "Friction") | (comparison["Test"] == "Reference")]
    velocity = comparison[(comparison["Test"] == "Velocity") | (comparison["Test"] == "Reference")]
    serpentinization = comparison[(comparison["Test"] == "Serpentinization") | (comparison["Test"] == "Reference")]

    names = ['viscosity', 'friction', 'velocity', 'serpentinization']

    for name in names:

        # Create the plot with 2 rows and 3 columns (for percentages and peak pressures)
        fig, ax = plt.subplots(2, 3, figsize=(15, 10))
        fig.suptitle('Exhumation and stagnation as a function of ' + name)

        # Plot percentages (First row)
        eval(name).plot(x='Model', y=['Tot_exh', 'Tot_stag'], kind='bar', ax=ax[0, 0], title='Percentage total exhumed particles')
        ax[0, 0].set_ylabel('Percentage')

        eval(name).plot(x='Model', y=['Sed_exh', 'Oc_exh', 'Ecl_exh'], kind='bar', ax=ax[0, 1], title='Percentage exhumed particles by lithology', colormap='tab20b')
        ax[0, 1].set_ylabel('Percentage')

        eval(name).plot(x='Model', y=['Sed_stag', 'Oc_stag', 'Ecl_stag'], kind='bar', ax=ax[0, 2], title='Percentage stagnant particles by lithology', colormap='tab20b')
        ax[0, 2].set_ylabel('Percentage')

        # Plot peak pressures (Second row) as line plots with markers
        models = eval(name)['Model'].unique()

        # Line plot for total exhumed and stagnant particles
        eval(name).plot(x='Model', y=['Exh_peakP', 'Stag_peakP'], kind='line', ax=ax[1, 0], title='Peak pressure total exhumed particles', marker='o')
        ax[1, 0].set_ylabel('Peak pressure (GPa)')
        ax[1, 0].set_xticks(range(len(models)))
        ax[1, 0].set_xticklabels(models)

        # Line plot for exhumed particles by lithology
        eval(name).plot(x='Model', y=['Exh_sed_peakP', 'Exh_oc_peakP', 'Exh_ecl_peakP'], kind='line', ax=ax[1, 1], title='Peak pressure exhumed particles by lithology', marker='o', colormap='tab20b')
        ax[1, 1].set_ylabel('Peak pressure (GPa)')
        ax[1, 1].set_xticks(range(len(models)))
        ax[1, 1].set_xticklabels(models)

        # Line plot for stagnant particles by lithology
        eval(name).plot(x='Model', y=['Stag_sed_peakP', 'Stag_oc_peakP', 'Stag_ecl_peakP'], kind='line', ax=ax[1, 2], title='Peak pressure stagnant particles by lithology', marker='o', colormap='tab20b')
        ax[1, 2].set_ylabel('Peak pressure (GPa)')
        ax[1, 2].set_xticks(range(len(models)))
        ax[1, 2].set_xticklabels(models)

       
        # Adjust layout and save the plot
        plt.tight_layout()
        plt.savefig("combined/peakP_" + name + ".png", dpi=500)
        plt.close()

    
if __name__ == '__main__':
    main()
