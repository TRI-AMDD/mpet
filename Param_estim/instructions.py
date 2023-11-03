# instruction.py contains the instructions for the parameter estimation
import numpy as np
import os

def sim_instructions():
    # define lower and upper bound for each parameter
    ensable_system = [
            [("Conductivity", "G_mean_c"), ["5e-16", "1e-14"]],
            [("Conductivity","sigma_s_c"), ["0.5","0.8"]],
            [("Electrodes", "k0_foil"), ["10","100"]],
            # [("Particles","mean_c"), ["90e-9","120e-9"]],
            # [("Particles","stddev_c"), ["20e-9","40e-9"]],
                    ]

    ensable_cathode = [
            [("Reactions", "k0"), ["0.5","2"]],
                    ]

    ensable_anode = []

    # Define C-rates and Temperatures
    c_rates = [1,2,5]  # C-rates

    temperatures = [298]  # K

    # Define params_system file
    params_system = 'params_system.cfg'
    # Define materials files
    param_c = 'params_LFP.cfg'
    param_a = 'params_graphite.cfg'  # only used if Nvol_a != 0, to be defined anyway

    parameter_ranges = [ensable_system, ensable_cathode, ensable_anode]
    operating_conditions = (c_rates, temperatures)
    config_files = [params_system, param_c, param_a]

    return parameter_ranges, operating_conditions, config_files


def norm_crate():
    """
    Normalization factor for the C-rates and capacity.
    """
    return 1

def weight_mse(temp, crate):
    """
    Assigne different weights depending on temperature and C-rate.
    """
    weight_mse = 1
    return weight_mse


def discharge():
    """
    Select if simulating charge or discharge.
    """
    discharge = True
    return discharge

def initial_c():
    """
    select initial concentration
    """
    c0 = 0
    return c0

def neglect_parts_volts():
    """
    datapoints of the voltage curve to neglect at the beginning of the simulation
    """
    nbr_points_init_negl = 1
    return nbr_points_init_negl


def data_path():
    """
    Experimental data path
    """
    data_path = r'C:\Users\file.csv'
    return data_path


def mpet_folder():
    """
    MPET folder
    """
    mpet_folder = r'C:\Users\mpet'
    return mpet_folder


def storage_folder():
    """
    Folder to store results into MPET folder
    """
    store = r'store'
    store_folder = os.path.join(mpet_folder(), store, 'prm_est_results')
    return store_folder
