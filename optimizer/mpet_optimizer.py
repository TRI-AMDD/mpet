#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import os
import shutil
import pickle
from scipy import integrate
from scipy import interpolate
from scipy import optimize
import subprocess
import sys

import mpet.outmat2txt as outmat2txt
import mpet.main as main
import mpet.io_utils as IO

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


def f(inputs, expt_QV, A, params_sys, params_cathode, params_anode, params_list, f_handle, x_0):
    """The function that finds the residual of the input parameter values.
       Takes in the values of the parameters that we are optimizing in inputs.
       args is [pickle_file_name, area of electrode(m^2), params_sys, params_cathode, params_anode,
       [list of params to optimize], file_handle of iteration file]
       replaces them, and runs
       an MPET simulation. Outputs the residual between experimental data and MPET simulation
       by calling function plot_area_difference"""
    #if sim_output already exists, remove
    if os.path.exists("sim_output") and os.path.isdir("sim_output"):
        shutil.rmtree("sim_output")
    #sed the parameters that we are trying to sub pythonically
    modify_mpet_params(params_sys, params_cathode, params_anode, params_list, np.multiply(inputs, x_0))
    main.main(params_sys) #runs simulation after subbing new parameters
    #run simulation, call plotter to get text file
    #if output_data.mat doesn't exist, then we set residual to large number because simulation failed
    residual = 0
    if os.path.isfile("sim_output/output_data.mat") and os.stat("sim_output/output_data.mat").st_size == 0:
        outmat2txt.main("sim_output")  #generates data
        general_data = np.loadtxt("sim_output/generalData.txt", delimiter = ",")
        #now saves data in general_data.txt
        mpet_Q, mpet_V = calculate_VQ(general_data, A) #generates MPET Q, V data
        residual = area_difference(mpet_Q, mpet_V, expt_QV[0,:], expt_QV[1,:]) #find residual of plot area
        output_string = "Iteration: " + str(inputs) + "  residual:  " + str(residual) + '\n'
        f_handle.write(output_string)
    else:
        residual = 1e8
        f_handle.write("Simulation failed with inputs:" + str(inputs) + '\n')
    print(residual)
    return residual


def modify_mpet_params(params_sys, params_cathode, params_anode, params_to_modify, param_values):
    """Modifies the inputs of the mpet parameter files.
       Inputs are the parameter file, a list of the parameters of modify,
       and respectively the parameter values to modify them to. Cannot
       put in dict structure because scipy optimize is weird.
       No output, just returns the modified files"""
    #use subprocess.run()
    for index in range(len(params_to_modify)):
        #sed in params_sys 
        param_name = params_to_modify[index]
        file_name = params_sys
        if params_to_modify[index][-7:] == "cathode":
             #if in cathode, we modify the corresponding params
             param_name = params_to_modify[index][:-8]
             file_name = params_cathode
             #sed -i "s/$i*/c\$i = $j" params_sys
        elif params_to_modify[index][-5:] == "anode":
             #if in anode, modify corresponding params
             param_name = params_to_modify[index][:-6]
             file_name = params_anode
        args = "sed -i '/" + param_name + " = "+ "*/c\\" + param_name + " = " + str(param_values[index]) + "' " + file_name
        subprocess.run([args], shell = True) #runs bash command through python
    return


def calculate_VQ(data, A):
    """Calculates the Q (in Ah), then multiplies by area to find the capacity.
       Inputs are generalData np array from the results of our simulation, and the area in m^2
       Output is the array of voltages and capacities during filling"""
    V_array = data[:,3]
    time_array = data[:,0]/3600 #converts s to hr
    I_array = data[:,5]
    #Q = int(i*dt) * A in Ah/m^2 * m^2
    Q_array = integrate.cumtrapz(I_array, time_array, initial  = 0) * A
    return Q_array, V_array


def area_difference(x1, y1, x2, y2):
    """Finds the difference between two curves (x1, y1) and (x2, y2).
        Returns the area difference. Is a metric to compare the two results, but can possibly
        switch metrics in the future."""
    fun1 = interpolate.interp1d(x1, y1, bounds_error = False, fill_value = 0)
    fun2 = interpolate.interp1d(x2, y2, bounds_error = False, fill_value = 0)
    max_x = max([np.amax(x1), np.amax(x2)])
    min_x = min([np.amin(x1), np.amin(x2)])
    xnew = np.linspace(min_x, max_x, 400)
    area_diff = np.trapz(np.abs(fun1(xnew)-fun2(xnew)), x = xnew)
    return area_diff


def process_experimental_data(pickle_file_name):
    """Reloads experimental data pickle file. inputs pickle file name for a batch
       returns the Q, V of the pickle file for a specific cell"""
    bat_dict = pickle.load(open(pickle_file_name, 'rb'))
    Q_array = bat_dict['b3c43']['cycles']['10']['Qd']
    V_array = bat_dict['b3c43']['cycles']['10']['V']
    QV_array = np.vstack((Q_array, V_array))
    return QV_array


def plot_functions(opt_args, expt_QV, A):
    """Plots the final result of the data and the optimal function.
       Inputs: optimal_arguments, data.
       No outputs."""
    #run simulation, call plotter to get text file
    outmat2txt.main("sim_output")  #generates data
    general_data = np.loadtxt("sim_output/generalData.txt", delimiter = ",")
    #now saves data in general_data.txt
    mpet_Q, mpet_V = calculate_VQ(general_data, A)
    plt.plot(mpet_Q, mpet_V, 'r', expt_QV[:,0], expt_QV[:,1], 'b')
    plt.xlabel('Q (A*hr)')
    plt.ylabel('V (V)')
    plt.legend(['MPET Simulation Fit', 'Experimental Data'])
    plt.savefig('final_fit.png')
    return


#puts in initial guesses and parameters to tweak
params_list = ['k0_cathode', 'k0_anode', 'alpha_cathode', 'alpha_anode', 'D_cathode', 'D_anode']
x0 = [1.6e-1, 3.0e+1, 0.5, 0.5, 5.3e-19, 1.25e-12]

#we feed in x0prime = all ones after rescaled
x0_prime = np.ones(6)

#process and extract experimental data
expt_QV_dat = process_experimental_data(pickle_name)

#tests to see if parameter file works!! 
try:
    main.main(mpet_params_sys)
except IndexError:
    print("ERROR: No parameter file specified. Aborting")
    raise

#open iteration file
f_handle = open('iteration_file', 'w')

#generates function arguments and variable bounds for parameters
arg_list = (expt_QV_dat, area, mpet_params_sys, mpet_params_cathode, mpet_params_anode, params_list, f_handle, x0)
b1 = np.zeros(len(x0))
b2 = np.array([np.inf, np.inf, 1/x0[2], 1/x0[3], np.inf, np.inf])
A = np.eye(len(x0))

cons = [{"type": "ineq", "fun": lambda x: A@x - b1}, {"type": "ineq", "fun": lambda x: -A@x + b2}]
result = optimize.minimize(f, x0_prime, args = arg_list, constraints = cons, method = 'COBYLA', options = {'rhobeg': 1})
print(result.success) # check if solver was successful

print("Final Optimal Parameters", np.multiply(result.x, x0))

plot_functions(np.multiply(result.x, x0), expt_QV_dat, area)
