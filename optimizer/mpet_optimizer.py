#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import os
import pickle
from scipy import integrate
from scipy import optimize
import subprocess
import sys

import mpet.outmat2txt as outmat2txt
import mpet.main as main

############Argument Parser & Docstrings#######################

if len(sys.argv) != 6:
    print("Usage: ./optimizer.py params_system.cfg params_cathode.cfg params_anode.cfg batch.pkl area_of_electrodes(m^2)")
    sys.exit()

mpet_params_sys = sys.argv[1]
mpet_params_cathode = sys.argv[2]
mpet_params_anode = sys.argv[3]
pickle_name = sys.argv[4]
area = float(sys.argv[5]) #in m^2

############Functions###########################

#parameters in the file that we choose to optimize:
#if they are in the sys file, just use the name of the parameter
#if they in the cathode or anode parameter file, add "_anode/_cathode"
#to the end of the parameter


def f(inputs, pickle_file, A, params_sys, params_cathode, params_anode, params_list):
    """The function that finds the residual of the input parameter values.
       Takes in the values of the parameters that we are optimizing in inputs.
       args is [pickle_file_name, area of electrode(m^2), params_sys, params_cathode, params_anode,
       [list of params to optimize]]
       replaces them, and runs
       an MPET simulation. Outputs the residual between experimental data and MPET simulation
       by calling function plot_area_difference"""
    #sed the parameters that we are trying to sub pythonically
    modify_mpet_params(params_sys, params_cathode, params_anode, params_list, inputs)
    main.main(params_sys)
    #run simulation, call plotter to get text file
    outmat2txt.main("sim_output") 
    general_data = np.loadtxt("sim_output/generalData.txt", delimiter = ",")
    #now saves data in general_data.txt
    mpet_Q, mpet_V = calculate_VQ(general_data, A)
    expt_Q, expt_V = process_experimental_data(pickle_file_name)
    residual = area_difference(mpet_Q, mpet_V, expt_Q, expt_V) #find residual of plot area
    return residual


def modify_mpet_params(param_sys, params_cathode, params_anode, params_to_modify, param_values):
    """Modifies the inputs of the mpet parameter files.
       Inputs are the parameter file, a list of the parameters of modify,
       and respectively the parameter values to modify them to. Cannot
       put in dict structure because scipy optimize is weird.
       No output, just returns the modified files"""
    #use subprocess.run()
    for index in range(len(params_to_modify)):
        args = "sed -i '/" + params_to_modify[index] + "*/c\\" + params_to_modify[index] + " = " + param_values[index] + " "
        if params_to_modify[index][-7:] == "cathode":
            args = args + params_cathode
            #sed -i "s/$i*/c\$i = $j" params_sys
        elif params_to_modify[index][-5:] == "anode":
            args = args + params_anode
        else:
            #sed in params_sys 
            args = args + params_sys
        subprocess.run(args)
    return


def calculate_VQ(data, A):
    """Calculates the Q (in Ah/m^2), then multiplies by area to find the capacity.
       Inputs are generalData np array from the results of our simulation, and the area in m^2
       Output is the array of voltages and capacities during filling"""
    V_array = data[:,3]
    time_array = data[:,0]
    I_array = data[:,5]
    Q_array = integrate.cumtrapz(I_array, time_array) * A
    return Q_array, V_array


def area_difference(x1, y1, x2, y2):
    """Finds the difference between two curves (x1, y1) and (x2, y2).
        Returns the area difference. Is a metric to compare the two results, but can possibly
        switch metrics in the future."""
    area_diff = np.abs(np.trapz(y2, x2) - np.trapz(y1, x1))
    return area_diff


def process_experimental_data(pickle_file_name):
    """Reloads experimental data pickle file. inputs pickle file name for a batch
       returns the Q, V of the pickle file for a specific cell"""
    bat_dict = pickle.load(open(pickle_name, 'rb'))
    Q_array = bat_dict['b3c43']['cycles']['10']['Qd']
    V_array = bat_dict['b3c43']['cycles']['10']['V']
    return Q_array, V_array


#tests to see if parameter file works!! 
try:
    main.main(mpet_params_sys)
except IndexError:
    print("ERROR: No parameter file specified. Aborting")
    raise

#puts in initial guesses and parameters to tweak
params_list = ['k0_cathode', 'k0_anode', 'alpha_cathode', 'alpha_anode', 'D_cathode', 'D_anode']
x0 = [1.6e-1, 3.0e-2, 0.5, 0.5, 5.3e-19, 1.25e-12]

arg_list = [pickle_name, area, mpet_params_sys, mpet_params_cathode, mpet_params_anode, params_list]

result = optimize.minimize(f, x0, args = arg_list, method = 'BFGS')
print(result.success) # check if solver was successful
