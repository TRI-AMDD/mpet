import subprocess
import shutil
import configparser
from instructions import *
import scipy.io as sio
import h5py
import numpy as np
import pickle
import copy
import pandas as pd
import random


def import_data(temps, protocols, c_rates_cc, c_rates_gitt, cont_volts):
    # csv file
    # first raw is header
    # header = "Temp(K)\tCondition(C-Rate/Volts)\tSpecCap(mAh/cm2)\tStepTime(s)\tOutput(V/C-rate)\n"

    dict_data = {}
    for protocol in protocols:
        dict_data[protocol] = {}
        data_path = data_files(protocol=protocol)
        mpet_folder_path = mpet_folder()
        path = os.path.join(mpet_folder_path,'Param_estim' ,data_path)
        
        data = pd.read_csv(path, header=0, sep="\t")
        # select data for the specified temperatures and C-rates
        # data = data.loc[data['Temp(K)'].isin(temps)]
        
        if protocol == 'cc':
            applied_conditions = c_rates_cc
        elif protocol == 'gitt':
            applied_conditions = c_rates_gitt
        elif protocol == 'pitt':
            applied_conditions = cont_volts

        # data = data.loc[data['Condition(C-Rate/Volts)'].isin(applied_conditions)]
        # convert data to dictionary
        t_ind = 0
        for temp in temps:
            dict_data[protocol][temp] = {}
            data_temp = data.loc[data['Temp(K)'] == temp]
            for cond in applied_conditions[t_ind]:
                dict_data[protocol][temp][cond] = {}
                data_crate = data_temp.loc[data_temp['Condition(C-Rate/Volts)'] == cond]
                dict_data[protocol][temp][cond]['Output(V/C-rate)'] = np.array(data_crate['Output(V/C-rate)'])
                dict_data[protocol][temp][cond]['StepTime(s)'] = np.array(data_crate['StepTime(s)'])
                dict_data[protocol][temp][cond]['SpecCap(mAh/cm2)'] = np.array(data_crate['SpecCap(mAh/cm2)'])
            t_ind += 1
    return dict_data



def call_process(mpet_dir, params_system):
    params_system_dir = os.path.join(mpet_dir,'configs', params_system)
    run_dir = os.path.join(mpet_dir,'bin', 'mpetrun.py')
    subprocess.call(["python", run_dir, params_system_dir])

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

def get_voltage(sim_output, file_type, mpet_dir):
    T_ref = 298.
    k = 1.381e-23
    e = 1.602e-19
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

def get_ff(sim_output, file_type, mpet_dir):
    if file_type == 'mat':
        ff = sim_output['ffrac_c'][0]
    elif file_type == 'hdf5':
        ff = np.squeeze(sim_output['ffrac_c'][...])
    return ff

def get_specap(sim_output, file_type, mpet_dir):
    if file_type == 'mat':
        ff = sim_output['ffrac_c'][0]
    elif file_type == 'hdf5':
        ff = np.squeeze(sim_output['ffrac_c'][...])
    specap = ff*spec_cap()
    return specap

def get_time(sim_output, file_type, mpet_dir):
    derived_values_pickle = os.path.join(mpet_dir,r'sim_output','input_dict_derived_values.p')
    with open(derived_values_pickle, 'rb') as f:
        dict_derived_values = pickle.load(f)
    td = dict_derived_values["t_ref"]
    if file_type == 'mat':
        time = sim_output['phi_applied_times'][0]*td 
    elif file_type == 'hdf5':
        time = np.squeeze(sim_output['time'][...])
    return time

def take_data_sim_out(mpet_dir):
    """
    Extract Voltage, Capacity and Times curves from the folders.
    """
    sim_output, file_type = get_sim_out_data(mpet_dir)
    volt = get_voltage(sim_output, file_type,mpet_dir)
    spec_cap = get_specap(sim_output, file_type,mpet_dir)
    time = get_time(sim_output, file_type,mpet_dir)
    return volt, spec_cap, time

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

def build_gitt_segments(c_rate_gitt, numb_time_steps, ocv_time):
    """
    Build the GITT segments.
    """

    time_per_step = 60/(c_rate_gitt*numb_time_steps)
    segments = "["
    for i in range(numb_time_steps-5):
        string = f"({c_rate_gitt},{time_per_step}), (0,{ocv_time}),"
        segments += string
    segments += "]"
    return segments


def update_params_system(ensambles, params_system, param_c, param_a, protocol, temp, condition):
    """
    Update the cfg files with the parameters in ensambles.
    """
    param_system_path = os.path.join(mpet_folder(),'configs', params_system)
    param_c_path = os.path.join(mpet_folder(),'configs', param_c)
    param_a_path = os.path.join(mpet_folder(),'configs', param_a)
    cfg_s = configparser.ConfigParser()
    cfg_s.optionxform = str
    cfg_s.read(param_system_path)
    new_params_system = cfg_s
    ensambles_s = ensambles[0]
    new_params_system["Electrodes"]["cathode"] = param_c
    new_params_system["Electrodes"]["anode"] = param_a
    if protocol == 'cc':
        new_params_system["Sim Params"]["profileType"] = 'CC'
        new_params_system["Sim Params"]["Crate"] = str(condition)
        new_params_system["Sim Params"]["T"] = str(temp)
        new_params_system["Sim Params"]["tramp"] = "0.01"
    elif protocol == 'gitt':
        new_params_system["Sim Params"]["profileType"] = 'CCsegments'
        new_params_system["Sim Params"]["segments"] = build_gitt_segments(condition, 10, 30)
        new_params_system["Sim Params"]["T"] = str(temp)
        new_params_system["Sim Params"]["tramp"] = "1"
    elif protocol == 'pitt':
        new_params_system["Sim Params"]["profileType"] = 'CV'
        new_params_system["Sim Params"]["V"] = str(condition)
        new_params_system["Sim Params"]["tend"] = "1e3" # minutes
        new_params_system["Sim Params"]["T"] = str(temp)
        new_params_system["Sim Params"]["tramp"] = "0.001"
    for val in ensambles_s:
        section = val[0][0]
        key = val[0][1]
        value = val[1][0]
        new_params_system[section][key] = value
        
        

    with open(param_system_path, 'w') as f:
        new_params_system.write(f)
    
    cfg_c = configparser.ConfigParser()
    cfg_c.optionxform = str
    cfg_c.read(param_c_path)
    new_param_c = cfg_c
    ensambles_c = ensambles[1]
    for val in ensambles_c:
        section = val[0][0]
        key = val[0][1]
        value = val[1][0]
        new_param_c[section][key] = value
    with open(param_c_path, 'w') as f:
        new_param_c.write(f)

    # if anode is present
    if len(ensambles[2]) != 0:
        cfg_a = configparser.ConfigParser()
        cfg_a.optionxform = str
        cfg_a.read(param_a_path)
        new_param_a = cfg_a
        ensambles_a = ensambles[2]
        for val in ensambles_a:
            section = val[0][0]
            key = val[0][1]
            value = val[1][0]
            new_param_a[section][key] = value
        with open(param_a_path, 'w') as f:
            new_param_a.write(f)

def take_original_tolerance(params_system):
    """
    Take the original tolerance from the cfg files.
    """
    param_system_path = os.path.join(mpet_folder(),'configs', params_system)
    cfg_s = configparser.ConfigParser()
    cfg_s.optionxform = str
    cfg_s.read(param_system_path)
    relTol = float(cfg_s.get('Sim Params', 'relTol'))
    absTol = float(cfg_s.get('Sim Params', 'absTol'))
    return relTol, absTol

def take_original_discretization(param_c):
    """
    Take the original discretization from the cfg files.
    """
    param_c_path = os.path.join(mpet_folder(),'configs', param_c)
    cfg_c = configparser.ConfigParser()
    cfg_c.optionxform = str
    cfg_c.read(param_c_path)
    discretization = float(cfg_c.get('Particles', 'discretization'))
    return discretization


def update_discretization(ensambles, param_c):
    """
    Update the cfg files with the lower discretization.
    """
    param_c_path = os.path.join(mpet_folder(),'configs', param_c)
    cfg_c = configparser.ConfigParser()
    cfg_c.optionxform = str
    cfg_c.read(param_c_path)
    new_params_system = cfg_c
    original_discretization = take_original_discretization(param_c)
    ensambles_c = ensambles[1]
    for val in ensambles_c:
        section = val[0][0]
        key = val[0][1]
        value = val[1][0]
        new_params_system[section][key] = value
        new_params_system["Particles"]["discretization"] = str(1e-9)
    with open(param_c_path, 'w') as f:
        new_params_system.write(f)


def reset_discretization(ensambles, param_c, discr):
    """
    Reset the cfg files with the original discretization.
    """
    param_c_path = os.path.join(mpet_folder(),'configs', param_c)
    cfg_c = configparser.ConfigParser()
    cfg_c.optionxform = str
    cfg_c.read(param_c_path)
    new_params_system = cfg_c
    ensambles_c = ensambles[1]
    for val in ensambles_c:
        section = val[0][0]
        key = val[0][1]
        value = val[1][0]
        new_params_system[section][key] = value
        new_params_system["Particles"]["discretization"] = str(discr)
    with open(param_c_path, 'w') as f:
        new_params_system.write(f)

    


def update_tolerance(ensambles, params_system):
    """
    Update the cfg files with the lower tolerance.
    """
    params_system_path = os.path.join(mpet_folder(),'configs', params_system)
    cfg_s = configparser.ConfigParser()
    cfg_s.optionxform = str
    cfg_s.read(params_system_path)
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
    with open(params_system_path, 'w') as f:
        new_params_system.write(f)

def reset_tolerance(ensambles, params_system, relTol, absTol):
    """
    Reset the cfg files with the original tolerance.
    """
    params_system_path = os.path.join(mpet_folder(),'configs', params_system)
    cfg_s = configparser.ConfigParser()
    cfg_s.optionxform = str
    cfg_s.read(params_system_path)
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
    with open(params_system_path, 'w') as f:
        new_params_system.write(f)




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





def mat_param_name(ensambles_tot_new):
    # create a name from the parameters
    nicename_mat = []
    for ensambles in ensambles_tot_new:
        # if anode is not present
        if len(ensambles) == 0:
            continue
        for val in ensambles:
            # convert val[1][0] in float, remove excess numbers after comma 
            # and convert in string
            float_val = float(val[1][0])
            float_val = "{:.2e}".format(float_val)
            value_name = str(float_val)
            nicename_mat.append(val[0][1] + "=" + str(value_name))
    nicename_mat = "_".join(nicename_mat)
    return nicename_mat

def make_folders(ensambles_tot_new, protocols, temperatures):
    # create store folder
    save_fold = save_folder()
    store_fold_mat = mat_param_name(ensambles_tot_new)
    store_fold_mat = os.path.join(save_fold, store_fold_mat)
    # create store folder for that parameter set
    if not os.path.exists(store_fold_mat):
        os.makedirs(store_fold_mat)
    for prot in protocols:
        for temps in temperatures:
            store_folde_temp = os.path.join(store_fold_mat, prot ,"T=" + str(temps))
            # create store folder for that temperature
            if not os.path.exists(store_folde_temp):
                os.makedirs(store_folde_temp)
    return store_fold_mat

def generate_tuple_bounds(ensamble):
    parameter_tuples = []
    for ensables in ensamble:
        for param_tuple in ensables:
            if len(param_tuple[1]) == 2:
                lower_limit, upper_limit = param_tuple[1]
                parameter_tuples.append((float(lower_limit), float(upper_limit)))
            else:
                raise Exception('Parameter range has to be of length 2')

    return tuple(parameter_tuples)


def generate_initial_guess(tuple_bounds, rand = True, seed = 1):
    random.seed(seed)
    initial_guess = []
    for bounds in tuple_bounds:
        if rand:
            initial_guess.append(random.uniform(bounds[0], bounds[1]))
        else:
            initial_guess.append((bounds[0] + bounds[1])/2)
    return initial_guess