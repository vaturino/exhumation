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



    #######################################################
    ###### Plot velocity comparison by Peak Pressure ######
    #######################################################
    fig, ax = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle('Peak pressure conditions as a function of velocity')

    # Plot for total exhumed and stagnant particles
    veloc.plot(x='Model', y=['Exh_peakP', 'Stag_peakP'], kind='line', ax=ax[0], title='Peak pressure', marker='o')
    ax[0].set_ylabel('Peak pressure (GPa)')
    ax[0].legend(loc='best')

    # Plot for exhumed particles by lithology
    veloc.plot(x='Model', y=['Exh_sed_peakP', 'Exh_oc_peakP', 'Exh_ecl_peakP'], kind='line', ax=ax[1], title='Peak pressure by lithology', marker='o', colormap='tab20b')
    ax[1].set_ylabel('Peak pressure (GPa)')
    ax[1].legend(loc='best')

    # Plot for stagnant particles by lithology
    veloc.plot(x='Model', y=['Stag_sed_peakP', 'Stag_oc_peakP', 'Stag_ecl_peakP'], kind='line', ax=ax[2], title='Peak pressure by lithology', marker='o', colormap='tab20b')
    ax[2].set_ylabel('Peak pressure (GPa)')
    ax[2].legend(loc='best')

    plt.tight_layout()
    plt.savefig("peakP/PeakP_velocity_test.eps", dpi=1000)
    plt.close()


    ########################################################
    ###### Plot viscosity comparison by Peak Pressure ######
    ########################################################
    fig, ax = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle('Peak pressure conditions as a function of viscosity')

    # Plot for total exhumed and stagnant particles
    visc.plot(x='Model', y=['Exh_peakP', 'Stag_peakP'], kind='line', ax=ax[0], title='Peak pressure', marker='o')
    ax[0].set_ylabel('Peak pressure (GPa)')

    # Plot for exhumed particles by lithology
    visc.plot(x='Model', y=['Exh_sed_peakP', 'Exh_oc_peakP', 'Exh_ecl_peakP'], kind='line', ax=ax[1], title='Peak pressure by lithology', marker='o', colormap='tab20b')
    ax[1].set_ylabel('Peak pressure (GPa)')
    ax[1].legend(loc='best')

    # Plot for stagnant particles by lithology
    visc.plot(x='Model', y=['Stag_sed_peakP', 'Stag_oc_peakP', 'Stag_ecl_peakP'], kind='line', ax=ax[2], title='Peak pressure by lithology', marker='o', colormap='tab20b')
    ax[2].set_ylabel('Peak pressure (GPa)')
    ax[2].legend(loc='best')

    plt.tight_layout()
    plt.savefig("peakP/PeakP_viscosity_test.eps", dpi=1000)
    plt.close()


    #######################################################
    ###### Plot friction comparison by Peak Pressure ######
    #######################################################
    fig, ax = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle('Peak pressure conditions as a function of friction')

    # Plot for total exhumed and stagnant particles
    fric.plot(x='Model', y=['Exh_peakP', 'Stag_peakP'], kind='line', ax=ax[0], title='Peak pressure', marker='o')
    ax[0].set_ylabel('Peak pressure (GPa)')
    ax[0].legend(loc='best')
    
    # Plot for exhumed particles by lithology
    fric.plot(x='Model', y=['Exh_sed_peakP', 'Exh_oc_peakP', 'Exh_ecl_peakP'], kind='line', ax=ax[1], title='Peak pressure by lithology', marker='o', colormap='tab20b')
    ax[1].set_ylabel('Peak pressure (GPa)')
    ax[1].legend(loc='best')

    # Plot for stagnant particles by lithology
    fric.plot(x='Model', y=['Stag_sed_peakP', 'Stag_oc_peakP', 'Stag_ecl_peakP'], kind='line', ax=ax[2], title='Peak pressure by lithology', marker='o', colormap='tab20b')
    ax[2].set_ylabel('Peak pressure (GPa)')
    ax[2].legend(loc='best')

    plt.tight_layout()
    plt.savefig("peakP/PeakP_friction_test.eps", dpi=1000)
    plt.close()


    ###############################################################
    ###### Plot serpentinization comparison by Peak Pressure ######
    ###############################################################
    fig, ax = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle('Peak pressure conditions as a function of serpentinization')

    # Define model names
    models = serp['Model'].unique()

    # Plot for total exhumed and stagnant particles
    serp.plot(x='Model', y=['Exh_peakP', 'Stag_peakP'], kind='line', ax=ax[0], title='Peak pressure', marker='o')
    ax[0].set_ylabel('Peak pressure (GPa)')
    ax[0].set_xticks(range(len(models)))
    ax[0].set_xticklabels(models)
    ax[0].legend(loc='best')

    # Plot for exhumed particles by lithology
    serp.plot(x='Model', y=['Exh_sed_peakP', 'Exh_oc_peakP', 'Exh_ecl_peakP', 'Exh_serp_peakP'], kind='line', ax=ax[1], title='Peak pressure by lithology', marker='o', colormap='tab20b')
    ax[1].set_ylabel('Peak pressure (GPa)')
    ax[1].set_xticks(range(len(models)))
    ax[1].set_xticklabels(models)
    ax[1].legend(loc='best')

    # Plot for stagnant particles by lithology
    serp.plot(x='Model', y=['Stag_sed_peakP', 'Stag_oc_peakP', 'Stag_ecl_peakP', 'Stag_serp_peakP'], kind='line', ax=ax[2], title='Peak pressure by lithology', marker='o', colormap='tab20b')
    ax[2].set_ylabel('Peak pressure (GPa)')
    ax[2].set_xticks(range(len(models)))
    ax[2].set_xticklabels(models)
    ax[2].legend(loc='best')

    plt.tight_layout()
    plt.savefig("peakP/PeakP_serpentinization_test.eps", dpi=1000)
    plt.close()


    




    # # Plot velocity comparison by percentage of exhumed and stagnant particles
    # fig, ax = plt.subplots(1, 3, figsize=(15, 5))
    # fig.suptitle('Peak pressure conditions as a function of velocity')

    # veloc.plot(x='Model', y=['Exh_peakP', 'Stag_peakP'], kind='bar', ax=ax[0], title='Percentage total exhumed particles')
    # ax[0].set_ylabel('Peak pressure (GPa)')

    # veloc.plot(x='Model', y=['Exh_sed_peakP', 'Exh_oc_peakP', 'Exh_ecl_peakP'], kind='bar', ax=ax[1], title='Percentage exhumed particles by lithology', colormap='tab20b')
    # ax[1].set_ylabel('Peak pressure (GPa)')

    # veloc.plot(x='Model', y=['Stag_sed_peakP', 'Stag_oc_peakP', 'Stag_ecl_peakP'], kind='bar', ax=ax[2], title='Percentage stagnant particles by lithology', colormap='tab20b')
    # ax[2].set_ylabel('Peak pressure (GPa)')

    # plt.tight_layout()
    # plt.savefig("peakP/PeakP_velocity_test.eps", dpi=1000)
    # plt.close()


    # # Plot viscosity comparison by percentage of exhumed and stagnant particles
    # fig, ax = plt.subplots(1, 3, figsize=(15, 5))
    # fig.suptitle('Peak pressure conditions as a function of viscosity')

    # visc.plot(x='Model', y=['Exh_peakP', 'Stag_peakP'], kind='bar', ax=ax[0], title='Percentage total exhumed particles')
    # ax[0].set_ylabel('Peak pressure (GPa)')

    # visc.plot(x='Model', y=['Exh_sed_peakP', 'Exh_oc_peakP', 'Exh_ecl_peakP'], kind='bar', ax=ax[1], title='Percentage exhumed particles by lithology', colormap='tab20b')
    # ax[1].set_ylabel('Peak pressure (GPa)')

    # visc.plot(x='Model', y=['Stag_sed_peakP', 'Stag_oc_peakP', 'Stag_ecl_peakP'], kind='bar', ax=ax[2], title='Percentage stagnant particles by lithology', colormap='tab20b')
    # ax[2].set_ylabel('Peak pressure (GPa)')

    # plt.tight_layout()
    # plt.savefig("peakP/PeakP_viscosity_test.eps", dpi=1000)
    # plt.close()



    # # Plot friction comparison by percentage of exhumed and stagnant particles
    # fig, ax = plt.subplots(1, 3, figsize=(15, 5))
    # fig.suptitle('Peak pressure conditions as a function of friction')

    # fric.plot(x='Model', y=['Exh_peakP', 'Stag_peakP'], kind='bar', ax=ax[0], title='Percentage total exhumed particles')
    # ax[0].set_ylabel('Peak pressure (GPa)')

    # fric.plot(x='Model', y=['Exh_sed_peakP', 'Exh_oc_peakP', 'Exh_ecl_peakP'], kind='bar', ax=ax[1], title='Percentage exhumed particles by lithology', colormap='tab20b')
    # ax[1].set_ylabel('Peak pressure (GPa)')

    # fric.plot(x='Model', y=['Stag_sed_peakP', 'Stag_oc_peakP', 'Stag_ecl_peakP'], kind='bar', ax=ax[2], title='Percentage stagnant particles by lithology', colormap='tab20b')
    # ax[2].set_ylabel('Peak pressure (GPa)')

    # plt.tight_layout()
    # plt.savefig("peakP/PeakP_friction_test.eps", dpi=1000)
    # plt.close()


    # # Plot serpentinization comparison by percentage of exhumed and stagnant particles
    # fig, ax = plt.subplots(1, 3, figsize=(15, 5))
    # fig.suptitle('Peak pressure conditions as a function of serpentinization')

    # serp.plot(x='Model', y=['Exh_peakP', 'Stag_peakP'], kind='bar', ax=ax[0], title='Percentage total exhumed particles')
    # ax[0].set_ylabel('Peak pressure (GPa)')

    # serp.plot(x='Model', y=['Exh_sed_peakP', 'Exh_oc_peakP', 'Exh_ecl_peakP', 'Exh_serp_peakP'], kind='bar', ax=ax[1], title='Percentage exhumed particles by lithology', colormap='tab20b')
    # ax[1].set_ylabel('Peak pressure (GPa)')

    # serp.plot(x='Model', y=['Stag_sed_peakP', 'Stag_oc_peakP', 'Stag_ecl_peakP', 'Stag_serp_peakP'], kind='bar', ax=ax[2], title='Percentage stagnant particles by lithology', colormap='tab20b')
    # ax[2].set_ylabel('Peak pressure (GPa)')

    # plt.tight_layout()
    # plt.savefig("peakP/PeakP_serpentinization_test.eps", dpi=1000)
    # plt.close()


if __name__ == '__main__':
    main()