# run_mpet.py

import os
import subprocess
import shutil
import configparser
from load_data import *
from instructions import *
from utils_prm_est import *
import scipy.io as sio
import h5py
import numpy as np
import pickle
import copy


def call_process(mpet_dir, params_system):
    cwd = os.getcwd()
    os.chdir(mpet_dir)
    params_system_dir = os.path.join('configs', params_system)
    run_dir = os.path.join('bin', 'mpetrun.py')
    subprocess.call(["python", run_dir, params_system_dir])
    os.chdir(cwd)

def copy_sim_out_in_folder(mpet_dir, folder_path):
    """
    Copy the sim_output folder in the store folder.
    """
    sim_out_folder = os.path.join(mpet_dir, r'sim_output')
    # copy sim_out_folder to folder_path
    # if folder_path exists, delete it
    if os.path.exists(folder_path):
        shutil.rmtree(folder_path)
    shutil.copytree(sim_out_folder, folder_path)

def get_sim_out_data(mpet_dir):
    # Implement this function
    data_dir = os.path.join(mpet_dir, r'sim_output')
    if os.path.exists(os.path.join(data_dir, 'output_data.mat')):
        sim_output = sio.loadmat(os.path.join(data_dir, 'output_data.mat'))
        file_type = 'mat'
    elif os.path.exists(os.path.join(data_dir, 'output_data.hdf5')):
        sim_output = h5py.File(os.path.join(data_dir, 'output_data.hdf5'), 'r')
        file_type = 'hdf5'
    return sim_output, file_type

def get_voltage(mpet_dir):
    T_ref = 298.
    k = 1.381e-23
    e = 1.602e-19
    sim_output, file_type = get_sim_out_data(mpet_dir)
    if file_type == 'mat':
        volt0 = -(k*T_ref/e) * sim_output['phi_cell'][0]
    elif file_type == 'hdf5':
        volt0 = -(k*T_ref/e) * np.squeeze(sim_output['phi_cell'][...])
    else:
        raise ValueError("File type not recognized")
    
    anode_pickle = os.path.join(mpet_dir,r'sim_output','input_dict_anode.p')
    derived_values_pickle = os.path.join(mpet_dir,r'sim_output','input_dict_derived_values.p')

    with open(derived_values_pickle, 'rb') as f:
        dict_derived_values = pickle.load(f)

    Etheta_c = dict_derived_values["c"]["phiRef"]
    if os.path.exists(anode_pickle):
        Ethet_a = dict_derived_values["a"]["phiRef"]
    else:
        Ethet_a = 0

    volt = -(k*T_ref/e)*(Etheta_c - Ethet_a) + volt0

    return volt

def get_ff(mpet_dir):
    sim_output, file_type = get_sim_out_data(mpet_dir)
    if file_type == 'mat':
        ff = sim_output['ffrac_c'][0]
    elif file_type == 'hdf5':
        ff = np.squeeze(sim_output['ffrac_c'][...])
    return ff

def take_data_sim_out(mpet_dir):
    """
    Extract Voltage and capacity curves from the folders.
    """
    volt = get_voltage(mpet_dir)
    ff = get_ff(mpet_dir)
    return volt, ff

def cut_off_voltage(mpet_dir):
    """
    Cut off voltage set in cfg
    """
    cfg = configparser.ConfigParser()
    cfg.optionxform = str
    cfg.read(os.path.join(mpet_dir,'sim_output','input_params_system.cfg'))
    if discharge():
        cutoff_voltage = float(cfg.get('Sim Params', 'Vmin'))
    else:
        cutoff_voltage = float(cfg.get('Sim Params', 'Vmax'))
    return cutoff_voltage

def update_params_system(ensambles, params_system, param_c, param_a, temp, crate):
    """
    Update the cfg files with the parameters in ensambles.
    """
    cwd = os.getcwd()
    os.chdir(mpet_folder())
    cfg_s = configparser.ConfigParser()
    cfg_s.optionxform = str
    cfg_s.read(os.path.join('configs', params_system))
    new_params_system = cfg_s
    ensambles_s = ensambles[0]
    for val in ensambles_s:
        section = val[0][0]
        key = val[0][1]
        value = val[1][0]
        new_params_system[section][key] = value
        new_params_system["Electrodes"]["cathode"] = param_c
        new_params_system["Electrodes"]["anode"] = param_a
        new_params_system["Sim Params"]["T"] = str(temp)
        new_params_system["Sim Params"]["Crate"] = str(crate)
    with open(os.path.join('configs', params_system), 'w') as f:
        new_params_system.write(f)
    
    cfg_c = configparser.ConfigParser()
    cfg_c.optionxform = str
    cfg_c.read(os.path.join('configs', param_c))
    new_param_c = cfg_c
    ensambles_c = ensambles[1]
    for val in ensambles_c:
        section = val[0][0]
        key = val[0][1]
        value = val[1][0]
        new_param_c[section][key] = value
    with open(os.path.join('configs', param_c), 'w') as f:
        new_param_c.write(f)

    # if anode is present
    if len(ensambles[2]) != 0:
        cfg_a = configparser.ConfigParser()
        cfg_a.optionxform = str
        cfg_a.read(os.path.join('configs', param_a))
        new_param_a = cfg_a
        ensambles_a = ensambles[2]
        for val in ensambles_a:
            section = val[0][0]
            key = val[0][1]
            value = val[1][0]
            new_param_a[section][key] = value
        with open(os.path.join('configs', param_a), 'w') as f:
            new_param_a.write(f)
    os.chdir(cwd)

def take_original_tolerance(params_system):
    """
    Take the original tolerance from the cfg files.
    """
    cwd = os.getcwd()
    os.chdir(mpet_folder())
    cfg_s = configparser.ConfigParser()
    cfg_s.optionxform = str
    cfg_s.read(os.path.join('configs', params_system))
    relTol = float(cfg_s.get('Sim Params', 'relTol'))
    absTol = float(cfg_s.get('Sim Params', 'absTol'))
    os.chdir(cwd)
    return relTol, absTol


def update_tolerance(ensambles, params_system):
    """
    Update the cfg files with the lower tolerance.
    """
    cwd = os.getcwd()
    os.chdir(mpet_folder())
    cfg_s = configparser.ConfigParser()
    cfg_s.optionxform = str
    cfg_s.read(os.path.join('configs', params_system))
    new_params_system = cfg_s
    relTol, absTol = take_original_tolerance(params_system)
    ensambles_s = ensambles[0]
    for val in ensambles_s:
        section = val[0][0]
        key = val[0][1]
        value = val[1][0]
        new_params_system[section][key] = value
        new_params_system["Sim Params"]["relTol"] = str(relTol*1e-3)
        new_params_system["Sim Params"]["absTol"] = str(absTol*1e-3)
        # new_params_system["Electrodes"]["cathode"] = param_c
    with open(os.path.join('configs', params_system), 'w') as f:
        new_params_system.write(f)
    
    os.chdir(cwd)

def reset_tolerance(ensambles, params_system, relTol, absTol):
    """
    Reset the cfg files with the original tolerance.
    """
    cwd = os.getcwd()
    os.chdir(mpet_folder())
    cfg_s = configparser.ConfigParser()
    cfg_s.optionxform = str
    cfg_s.read(os.path.join('configs', params_system))
    new_params_system = cfg_s
    ensambles_s = ensambles[0]
    for val in ensambles_s:
        section = val[0][0]
        key = val[0][1]
        value = val[1][0]
        new_params_system[section][key] = value
        new_params_system["Sim Params"]["relTol"] = str(relTol)
        new_params_system["Sim Params"]["absTol"] = str(absTol)
        # new_params_system["Electrodes"]["cathode"] = param_c
    with open(os.path.join('configs', params_system), 'w') as f:
        new_params_system.write(f)
    
    os.chdir(cwd)



def scale(ensamble_tot):
    """
    Scale the parameters in ensamble
    """
    scaled_ensamble = []
    scaled_ensamble = copy.deepcopy(ensamble_tot)
    # scale the parameters
    scaling_vec = []
    for ensamble in scaled_ensamble:
        for param_tuple in ensamble:
            if len(param_tuple[1]) == 2:
                lower_limit, upper_limit = param_tuple[1]
                scaling = (float(upper_limit) - float(lower_limit))
                param_tuple[1] = [str(float(lower_limit)/scaling), str(float(upper_limit)/scaling)]
                scaling_vec.append(scaling)
            else:
                raise Exception('Parameter range has to be of length 2')
    return scaled_ensamble, scaling_vec

def rescale(params, scaling_vec):
    """
    Rescale the parameters in ensamble
    """
    parameters = []
    for value, scaling in zip(params, scaling_vec):
        parameters.append(value*scaling)
    return parameters
