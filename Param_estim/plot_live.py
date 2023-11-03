import os
import matplotlib.pyplot as plt
# from multiprocessing import Process
from instructions import *
import scipy.io as sio
import h5py
import numpy as np
import pickle
from load_data import *
import time


# Define a function to extract voltage and filling fraction from a folder (replace with your actual data extraction logic)
# Function to extract data from a folder

def get_sim_data(folder_path):
    # Implement this function
    data_dir = folder_path
    if os.path.exists(os.path.join(data_dir, 'output_data.mat')):
        sim_output = sio.loadmat(os.path.join(data_dir, 'output_data.mat'))
        file_type = 'mat'
    elif os.path.exists(os.path.join(data_dir, 'output_data.hdf5')):
        sim_output = h5py.File(os.path.join(data_dir, 'output_data.hdf5'), 'r')
        file_type = 'hdf5'
    return sim_output, file_type


def get_voltage_plot(folder_path):
    T_ref = 298.
    k = 1.381e-23
    e = 1.602e-19
    sim_output, file_type = get_sim_data(folder_path)
    if file_type == 'mat':
        volt0 = -(k*T_ref/e) * sim_output['phi_cell'][0]
    elif file_type == 'hdf5':
        volt0 = -(k*T_ref/e) * np.squeeze(sim_output['phi_cell'][...])
    else:
        raise ValueError("File type not recognized")
    
    anode_pickle = os.path.join(folder_path,'input_dict_anode.p')
    derived_values_pickle = os.path.join(folder_path,'input_dict_derived_values.p')

    with open(derived_values_pickle, 'rb') as f:
        dict_derived_values = pickle.load(f)

    Etheta_c = dict_derived_values["c"]["phiRef"]
    if os.path.exists(anode_pickle):
        Ethet_a = dict_derived_values["a"]["phiRef"]
    else:
        Ethet_a = 0

    volt = -(k*T_ref/e)*(Etheta_c - Ethet_a) + volt0

    return volt


def get_ff_plot(folder_path):
    sim_output, file_type = get_sim_data(folder_path)
    if file_type == 'mat':
        ff = sim_output['ffrac_c'][0]
    elif file_type == 'hdf5':
        ff = np.squeeze(sim_output['ffrac_c'][...])
    # normalize ff
    c0 = initial_c()
    # c0 = 0.88
    ff = (ff - c0)/(norm_crate())
    return ff


def extract_data(folder_path):
    dict_sim = {}
    folder_path = os.path.join(storage_folder(), folder_path)

    ensambles_tot, operating_conditions, config_files = sim_instructions()

    crates = operating_conditions[0]
    temperatures = operating_conditions[1]
    
    for temp in temperatures:
        dict_sim[temp] = {}
        path_temp = os.path.join(folder_path, "T=" + str(temp))
        for c in crates:
            path = os.path.join(path_temp, "C=" + str(c))
            dict_sim[temp][c] = {}
            dict_sim[temp][c]['V'] = get_voltage_plot(path)
            dict_sim[temp][c]['ff'] = get_ff_plot(path)

    return dict_sim


# Function to get the oldest and most recent folders within a directory
def get_oldest_and_recent_folders(directory_path):
    subfolders = [f for f in os.listdir(directory_path) if os.path.isdir(os.path.join(directory_path, f))]
    
    # Sort the subfolders by their creation time (oldest first)
    oldest_folder = min(subfolders, key=lambda x: os.path.getctime(os.path.join(directory_path, x)))

    # Sort the subfolders by their creation time (most recent first)
    # recent_folder = max(subfolders, key=lambda x: os.path.getctime(os.path.join(directory_path, x)))
    # second to last folder
    recent_folder = sorted(subfolders, key=lambda x: os.path.getctime(os.path.join(directory_path, x)))[-2]

    return oldest_folder, recent_folder


# Define a function to update the MSE vs. iterations plot
def update_mse_plot(mse_values, ax):
    ax.plot(mse_values, 'o', color = 'C0')
    ax.set_xlabel('Iterations')
    ax.set_ylabel('MSE')
    # ax.set_ylim([0, 0.5])


# Define a function to update the voltage curves plot
def update_voltage_curves(dict_old, dict_new, ax):
    ensambles_tot, operating_conditions, config_files = sim_instructions()
    crates = operating_conditions[0]
    temperatures = operating_conditions[1]
    dict_of_experiments = import_data(temperatures,crates)
    # plt.figure()
    # Update the voltage curves plot
    for temp in temperatures:
        for c in crates:
            voltage_old = dict_old[temp][c]['V']
            ff_old = dict_old[temp][c]['ff']
            voltage_new = dict_new[temp][c]['V']
            ff_new = dict_new[temp][c]['ff']

            vexp = dict_of_experiments[temp][c]['V']
            cexp = dict_of_experiments[temp][c]['NormCap']
            labels = 'T = ' + str(temp) + ', C = ' + str(c) + ', Exp'
            ax.plot(cexp, vexp, label=labels, marker='o', linestyle='-', color = 'black', markersize=2)

            # Plot the voltage curve of the first simulation
            labels = 'T = ' + str(temp) + ', C = ' + str(c) + ', First Sim'
            ax.plot(ff_old, voltage_old, label=labels, color = 'gray', linestyle='--')
            
            # Plot the voltage curve of the last simulation
            labels = 'T = ' + str(temp) + ', C = ' + str(c) + ', Last Sim'
            ax.plot(ff_new, voltage_new, label=labels)
    
    ax.set_xlabel('filling fraction')
    ax.set_ylabel('V')
    # ax.set_ylim([3.2, 4])
    ax.legend(fontsize=2)
    
    plt.pause(1)


# Define a function to run the simulation and update plots
def plt_simulation(): # Replace with the actual log file path
    folder_path = storage_folder()  # Replace with the actual folder path
    log_file = os.path.join(folder_path, 'log_book.txt')

    
    fig, ax = plt.subplots(1,2, figsize=(10, 5))
    num_of_lines = 0
    while True:
        # Simulate and log the results
        num_lines_old = num_of_lines
        # wait 10 seconds
        # time.sleep(10)
        mse, num_of_lines = read_mse(log_file)

        # Update the MSE vs. iterations plot
        # update only when a new line in log_file is added
        if num_of_lines > 2 and num_of_lines != num_lines_old:
            
            update_mse_plot(mse, ax[0])

            # Extract data from the latest simulation folder
            oldest_folder, recent_folder = get_oldest_and_recent_folders(folder_path)
            dict_old = extract_data(oldest_folder)
            dict_new = extract_data(recent_folder)

            update_voltage_curves(dict_old, dict_new, ax[1])

        # if log_book.txt contains "end" then stop
        with open(log_file, 'r') as f:
            lines = f.readlines()
            if lines[-1] == "end\n":
                # wait 30 seconds
                time.sleep(30)
                break
            else:
                plt.cla()
                continue


def read_mse(log_file):
    with open(log_file, 'r') as f:
        lines = f.readlines()
        mse_vec = []
        num_of_lines = 0
        for line in lines[1:]:
            if line != "final parameters\n":
                mse_value = line.split('\t')[-1]
                mse_value = float(mse_value.split('\n')[0])
                mse_vec.append(mse_value)
                num_of_lines += 1
            else:
                break
    return mse_vec, num_of_lines


# Create a process for running the simulation and updating plots
plt_simulation()

# Display the plots
plt.ion()  # Turn on interactive mode
plt.show()
