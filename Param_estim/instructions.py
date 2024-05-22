# instruction.py contains the instructions for the parameter estimation
import datetime
import os

def sim_instructions():
    # A lot of smalle 40nm particles with low wiring does not work 
    optimization_method = 'GD'  # 'DE' or 'GD'

    ensable_system = [
            [("Conductivity", "sig_carb_c"), ["1e-3", "5e-4"]],
            [("Conductivity", "E_sig_carb_c"), ["0.1", "0.4"]],
            # [("Electrodes", "k0_foil"), ["10","50"]],
                    ]

    ensable_cathode = [
            [("Reactions", "k0"), ["5", "30"]],
                    ]

    ensable_anode = []

    # Define C-rates and Temperatures
    
    # temps = [298,283,268]  # K
    temperatures = [298,283,268]  # K
    # Possible protocols: 'cc', 'pitt', 'gitt'
    protocols = ['cc']
    # one list for each protocol and temperature
    # it is important that the order of temperatures and conditions is the same
    c_rates_cc = [[1],
                  [2],
                  [0.5]]  # C-rates
    c_rates_gitt = [[1],
                    [1]]  # C-rates
    cont_volts = [[3.377, 3.351, 3.326],
                  [3.376,3.352,3.271,3.221]]  # V
    # volts at 298 = [3.376, 3.351 ,3.326 ,2.926 ,3.126 ,3.026]
    # volts at 268 = [3.376,3.352,3.271,3.221]

    # Define params_system file
    params_system = 'params_system_wiring.cfg'
    # Define param_parameters file
    param_c = 'params_LFP_wiring.cfg'
    param_a = 'params_graphite.cfg'  # only used if Nvol_a != 0

    parameter_ranges = [ensable_system, ensable_cathode, ensable_anode]
    operating_conditions = (protocols, c_rates_cc, c_rates_gitt, cont_volts, temperatures)
    config_files = [params_system, param_c, param_a]

    return parameter_ranges, operating_conditions, config_files, optimization_method

# assign different weights
def weight_mse(protocol, temp, applied_con):
    # function that can be personalized to 
    # assign different weights to different
    # operating conditions

    if protocol == 'pitt':
        weight_mse = 0.2
        # since the current goes from 0 to 10 and the voltage from 3.5 to 2.5
        # the error on the current must be reduced
    return weight_mse

def c0():
    c0 = 0.01
    return c0

def discharge():
    discharge = True
    return discharge

def sp_cap():
    spec_cap = 1 # mAh/cm2
    return spec_cap

def gitt_protocol():
    steps = 20
    rest_time = 30 # min
    # NB: the simulation will run for 70% of the real steps to avoid crashing
    return steps, rest_time
# if it is necessary to neglect the first
# initial points of the voltage curve
def neglect_parts_volts():
    nbr_points_init_negl = 1
    return nbr_points_init_negl

def data_files(protocol):
    data_path = r'exp_data'
    if protocol == 'cc':
        data_cc = os.path.join(data_path, 'cc.csv')
        return data_cc
    elif protocol == 'pitt':
        data_pitt = os.path.join(data_path, 'pitt.csv')
        return data_pitt
    elif protocol == 'gitt':
        data_gitt = os.path.join(data_path, 'gitt.csv')
        return data_gitt
    else:
        print('Protocol not recognized')
    return data_cc, data_pitt, data_gitt

def mpet_folder():
    mpet_folder = r'C:\Users\pierfrancescoo\Documents\Phase-field\mpet-LFMP\mpet'
    return mpet_folder

# folder = datetime.datetime.now().strftime("%Y-%m-%d_%H%M") + "_DE"

def save_folder():
    folder = os.path.join(mpet_folder(), 'Param_estim', "Results")
    return folder
