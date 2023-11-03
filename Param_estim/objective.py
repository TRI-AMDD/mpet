# mse.py
from util_run import *
from load_data import *
from utils_prm_est import *
from instructions import *
import numpy as np
from scipy.interpolate import interp1d
import time
import os
import copy


def run_mpet(scaled_parameters):
    """
    Create the folder structure to store the data from the simulations.
    Update the cfg file of params_system, params_c
    and params_a with parameters.
    Run the simulation.
    Copy the sim_output folder in the store folder.
    Extract Voltage and capacity curves from the folders.
    """
    
    ensamble_tot, operating_conditions, config_files = sim_instructions()
    scaled_enasmble, scaling_factors = scale(ensamble_tot)
    parameters = rescale(scaled_parameters, scaling_factors)

    crates = operating_conditions[0]
    temperatures = operating_conditions[1]

    ensambles_tot_new = []
    # use paramters and ensambles_tot to update the cfg files
    value_index = 0
    for parameter_list in ensamble_tot:
        for param_tuple in parameter_list:
            param_tuple[1] = [str(parameters[value_index])]
            value_index += 1
        ensambles_tot_new.append(parameter_list)
    store_fold_mat = make_folders(ensambles_tot_new, temperatures)

    dic_resul = {}
    datetime = os.path.getctime(os.path.join(mpet_folder(), 'sim_output'))
    for temps in temperatures:
        dic_resul[temps] = {}
        for c in crates:
            dic_resul[temps][c] = {}
            print('C rate = ', c)
            print('T = ', temps)
            store_folder_c = os.path.join(store_fold_mat, "T=" + str(temps), "C=" + str(c))

            norm_c = c*norm_crate()
            update_params_system(ensambles_tot_new,
                                 params_system = config_files[0],
                                 param_c = config_files[1],
                                 param_a = config_files[2],
                                 temp = temps,
                                 crate = norm_c)
            
            # wait 0.5 seconds to avoid overwriting
            time.sleep(0.5)
            call_process(mpet_folder(), config_files[0])

            creation_time = os.path.getctime(os.path.join(mpet_folder(), 'sim_output'))
            

            if creation_time == datetime:
                print('tol to high, simulation did not run')
                print('decresing tolerance and running again')

                relTol, absTol = take_original_tolerance(config_files[0])
                update_tolerance(ensambles_tot_new, config_files[0])
                call_process(mpet_folder(), config_files[0])
                creation_time = os.path.getctime(os.path.join(mpet_folder(), 'sim_output'))
                reset_tolerance(ensambles_tot_new, config_files[0], relTol, absTol)
                if creation_time == datetime:
                    raise ValueError("Simulation did not run")
                else:
                    volt, ff = take_data_sim_out(mpet_folder())
                    
            else:
                datetime = creation_time
                # take data
                volt, ff = take_data_sim_out(mpet_folder())
                
                
                # if simulation did not finish do it again with lower tolernace
                if ((discharge()==True
                    and (cut_off_voltage(mpet_folder()) < (volt[-1] - 0.2))) 
                    or ((discharge()==False)
                    and (cut_off_voltage(mpet_folder()) > (volt[-1] + 0.2)))):
                        print('simulation did not finish')
                        print('decresing tolerance and running again')
                        relTol, absTol = take_original_tolerance(config_files[0])
                        update_tolerance(ensambles_tot_new, config_files[0])
                        call_process(mpet_folder(), config_files[0])
                        reset_tolerance(ensambles_tot_new, config_files[0], relTol, absTol)
                        volt, ff = take_data_sim_out(mpet_folder())
                        if ((discharge()==True
                            and (cut_off_voltage(mpet_folder()) < (volt[-1] - 0.1))) 
                            or ((discharge()==False)
                            and (cut_off_voltage(mpet_folder()) > (volt[-1] + 0.1)))):
                            raise ValueError("Simulation did not finish")

            dic_resul[temps][c]["volt"] = volt
            dic_resul[temps][c]["ff"] = ff
            copy_sim_out_in_folder(mpet_folder(), store_folder_c)

    return dic_resul


def obj(parameters):
    dict_of_results = run_mpet(parameters)
    # remember you need to interpolate to have same points
    # interpolate simulation data to be able to compare with experimental data
    ensamble_tot, operating_conditions, config_files = sim_instructions()
    crates = operating_conditions[0]
    temperatures = operating_conditions[1]
    dict_of_experiments = import_data(temperatures,crates)
    weighted_mse = 0
    for temps in temperatures:
        for crate in crates:
            volt = dict_of_results[temps][crate]['volt']
            ff = dict_of_results[temps][crate]['ff']
            # normalize ff
            c0 = initial_c()
            ff = (ff - c0)/(norm_crate())

            vexp = dict_of_experiments[temps][abs(crate)]['V']
            cexp = dict_of_experiments[temps][abs(crate)]['NormCap']

            vinterp = np.interp(cexp, ff, volt)
            weighted_mse += np.average((vexp - vinterp)**2) * weight_mse(temps, crate)

    weighted_mse /= len(temperatures) * len(crates)
    print('round done')
    with open(os.path.join(storage_folder(),"log_book.txt"), 'a') as f:
        param_rescaled = copy.deepcopy(parameters)
        scaled_enasmble, scaling_factors = scale(ensamble_tot)
        params = rescale(param_rescaled, scaling_factors)
        for par in params:
            f.write(str(par) + "\t")
        f.write(str(weighted_mse) + "\n")

    return weighted_mse