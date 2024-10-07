#! /usr/bin/python3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def main():
    comparison = pd.read_excel('comparison.xlsx')

    visc = comparison[(comparison["Test"] == "Viscosity") | (comparison["Test"] == 'Reference')]
    fric = comparison[(comparison["Test"] == "Friction") | (comparison["Test"] == 'Reference')]
    veloc = comparison[(comparison["Test"] == "Velocity") | (comparison["Test"] == 'Reference')]
    serp = comparison[(comparison["Test"] == "Serpentinization") | (comparison["Test"] == 'Reference')]



    # Plot velocity comparison by percentage of exhumed and stagnant particles
    fig, ax = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle('Percentage of exhumed and stagnant particles as a function of velocity')

    veloc.plot(x='Model', y=['Tot_exh', 'Tot_stag'], kind='bar', ax=ax[0], title='Percentage total exhumed particles')
    ax[0].set_ylabel('Percentage')

    veloc.plot(x='Model', y=['Sed_exh', 'Oc_exh', 'Ecl_exh'], kind='bar', ax=ax[1], title='Percentage exhumed particles by lithology', colormap='tab20b')
    ax[1].set_ylabel('Percentage')

    veloc.plot(x='Model', y=['Sed_stag', 'Oc_stag', 'Ecl_stag'], kind='bar', ax=ax[2], title='Percentage stagnant particles by lithology', colormap='tab20b')
    ax[2].set_ylabel('Percentage')

    plt.tight_layout()
    plt.savefig("percentage/Percentage_velocity_test.eps", dpi=1000)
    plt.close()



    # Plot friction comparison by percentage of exhumed and stagnant particles
    fig, ax = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle('Percentage of exhumed and stagnant particles as a function of friction')

    fric.plot(x='Model', y=['Tot_exh', 'Tot_stag'], kind='bar', ax=ax[0], title='Percentage total exhumed particles')
    ax[0].set_ylabel('%')

    fric.plot(x='Model', y=['Sed_exh', 'Oc_exh', 'Ecl_exh'], kind='bar', ax=ax[1], title='Percentage exhumed particles by lithology', colormap='tab20b')
    ax[1].set_ylabel('%')

    fric.plot(x='Model', y=['Sed_stag', 'Oc_stag', 'Ecl_stag'], kind='bar', ax=ax[2], title='Percentage stagnant particles by lithology', colormap='tab20b')
    ax[2].set_ylabel('%')

    plt.tight_layout()
    plt.savefig("percentage/Percentage_friction_test.eps", dpi=1000)
    plt.close()


    # Plot viscosity comparison by percentage of exhumed and stagnant particles
    fig, ax = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle('Percentage of exhumed and stagnant particles as a function of viscosity')

    visc.plot(x='Model', y=['Tot_exh', 'Tot_stag'], kind='bar', ax=ax[0], title='Percentage total exhumed particles')
    ax[0].set_ylabel('%')

    visc.plot(x='Model', y=['Sed_exh', 'Oc_exh', 'Ecl_exh'], kind='bar', ax=ax[1], title='Percentage exhumed particles by lithology', colormap='tab20b')
    ax[1].set_ylabel('%')

    visc.plot(x='Model', y=['Sed_stag', 'Oc_stag', 'Ecl_stag'], kind='bar', ax=ax[2], title='Percentage stagnant particles by lithology', colormap='tab20b')
    ax[2].set_ylabel('%')

    plt.tight_layout()
    plt.savefig("percentage/Percentage_viscosity_test.eps", dpi=1000)
    plt.close()


    # Plot serpentinization comparison by percentage of exhumed and stagnant particles
    fig, ax = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle('Percentage of exhumed and stagnant particles as a function of serpentinization')

    serp.plot(x='Model', y=['Tot_exh', 'Tot_stag'], kind='bar', ax=ax[0], title='Percentage total exhumed particles')
    ax[0].set_ylabel('%')

    serp.plot(x='Model', y=['Sed_exh', 'Oc_exh', 'Ecl_exh', 'Serp_exh'], kind='bar', ax=ax[1], title='Percentage exhumed particles by lithology', colormap='tab20b')
    ax[1].set_ylabel('%')

    serp.plot(x='Model', y=['Sed_stag', 'Oc_stag', 'Ecl_stag', 'Serp_exh'], kind='bar', ax=ax[2], title='Percentage stagnant particles by lithology', colormap='tab20b')
    ax[2].set_ylabel('%')

    plt.tight_layout()
    plt.savefig("percentage/Percentage_serpentinization_test.eps", dpi=1000)
    plt.close()
    


if __name__ == '__main__':
    main()