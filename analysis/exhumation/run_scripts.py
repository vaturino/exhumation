#!/usr/bin/env python3

import subprocess

#List of model names
model_list = ['sed_bas', 'sed_bas_cttV', 'sed_bas_Vto0', 'sed_bas_Vto10pc', 'sed_bas_Vto50pc', 'sed_bas_halfd', 'sed_bas_doubled', 'sed_bas_mu0.2', 'sed_bas_mu0.05', 'sed_bas_diffFric', 'sed_bas_serp1km', 'sed_bas_serp_patch']

# List of programs to run
program_list = ['exhumation_timescales.py', 'stagnation_timescales.py', 'time_exh_stag.py']

# Loop through each program and call it with the JSON file
for model in model_list:
    print("Running model: " + model)
    for program in program_list:
        subprocess.call(['python3', program, model + '.json'])
        print("Finished: " + program)
