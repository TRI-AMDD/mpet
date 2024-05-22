from utils import *
from instructions import *
import numpy as np
import os
import copy
import time
import matplotlib.pyplot as plt

def run_mpet(scaled_parameters):
    """
    Create the folder structure to store the data from the simulations.
    Update the cfg file of params_system, params_c
    and params_a with parameters.
    Run the simulation.
    Copy the sim_output folder in the store folder.
    Extract Voltage and capacity curves from the folders.
    """
    
    ensamble_tot, operating_conditions, config_files, optimization_method = sim_instructions()
    # (protocols, c_rates_cc, c_rates_gitt, cont_volt, temperatures) = operating_conditions
    scaled_enasmble, scaling_factors = scale(ensamble_tot)
    parameters = rescale(scaled_parameters, scaling_factors)

    protocols = operating_conditions[0]
    c_rates_cc = operating_conditions[1]
    c_rates_gitt = operating_conditions[2]
    cont_volts = operating_conditions[3]
    temps = operating_conditions[4]

    ensambles_tot_new = []
    # use paramters and ensambles_tot to update the cfg files
    value_index = 0
    for parameter_list in ensamble_tot:
        for param_tuple in parameter_list:
            param_tuple[1] = [str(parameters[value_index])]
            value_index += 1
        ensambles_tot_new.append(parameter_list)
    store_fold_mat = make_folders(ensambles_tot_new, protocols, temps)

    dict_res = {}
    for prot in protocols:
        dict_res[prot] = {}
        for i, temp in enumerate(temps):
            dict_res[prot][temp] = {}
            if prot == 'cc':
                applied_conditions = c_rates_cc[i]
                fold_name = "C"
            elif prot == 'gitt':
                applied_conditions = c_rates_gitt[i]
                fold_name = "C"
            elif prot == 'pitt':
                applied_conditions = cont_volts[i]
                fold_name = "V"
            else:
                print('protocol not recognized')

            for cond in applied_conditions:
                print('running protocol: ', prot, ' at temperature: ', temp, ' at condition: ', cond)
                dict_res[prot][temp][cond] = {}
                store_folder = os.path.join(store_fold_mat, prot, "T=" + str(temp), str(cond)+fold_name)
                
                update_params_system(ensambles_tot_new,
                                    params_system = config_files[0],
                                    param_c = config_files[1],
                                    param_a = config_files[2],
                                    protocol = prot, 
                                    temp = temp, 
                                    condition = cond)
                
                # wait 1 second to avoid overwriting
                time.sleep(1)
                datetime = os.path.getctime(os.path.join(mpet_folder(), 'sim_output'))
                call_process(mpet_folder(), config_files[0])
                creation_time = os.path.getctime(os.path.join(mpet_folder(), 'sim_output'))
                

                if creation_time == datetime:
                    print('simulation did not run')                    
                    print('decreasing tolerance and running again')
                    relTol, absTol = take_original_tolerance(config_files[0])
                    update_tolerance(ensambles_tot_new, config_files[0])
                    call_process(mpet_folder(), config_files[0])
                    creation_time = os.path.getctime(os.path.join(mpet_folder(), 'sim_output'))
                    reset_tolerance(ensambles_tot_new, config_files[0], relTol, absTol)

                    if creation_time == datetime:
                        print('!!!!!! simulation did not run again with lower tolerance !!!!!!')
                        print('returning 0.25')
                        return 0.25
                    else:
                        volt, curr, spec_cap, tim = take_data_sim_out(mpet_folder())
                        
                else:
                    datetime = creation_time
                    # take data
                    volt, curr, spec_cap, tim = take_data_sim_out(mpet_folder())
                    
                    
                    # if simulation did not finish do it again with lower tolernace
                    if prot == 'cc':
                        voltage_tolarenace = 0.1
                        if ((discharge()==True
                            and (cut_off_voltage(mpet_folder()) < (volt[-1] - voltage_tolarenace))) 
                            or ((discharge()==False)
                            and (cut_off_voltage(mpet_folder()) > (volt[-1] + voltage_tolarenace)))):
                                print('simulation did not finish')
                                print('decresing tolerance and running again')

                                relTol, absTol = take_original_tolerance(config_files[0])
                                update_tolerance(ensambles_tot_new, config_files[0])
                                call_process(mpet_folder(), config_files[0])
                                reset_tolerance(ensambles_tot_new, config_files[0], relTol, absTol)
                                volt, curr, spec_cap, tim = take_data_sim_out(mpet_folder())
                                voltage_tolarenace*=1.5
                                if ((discharge()==True
                                    and (cut_off_voltage(mpet_folder()) < (volt[-1] - voltage_tolarenace))) 
                                    or ((discharge()==False)
                                    and (cut_off_voltage(mpet_folder()) > (volt[-1] + voltage_tolarenace)))):
                                    print('simulation did not finish again, returning 0.25')
                                    return 0.25
                    elif prot == 'gitt':
                        if spec_cap[-1] < 0.5*sp_cap():
                            print('simulation did not finish')
                            print('decreasing  tolerance and running again')
                            relTol, absTol = take_original_tolerance(config_files[0])
                            update_tolerance(ensambles_tot_new, config_files[0])
                            call_process(mpet_folder(), config_files[0])
                            reset_tolerance(ensambles_tot_new, config_files[0], relTol, absTol)
                            volt, curr, spec_cap, tim = take_data_sim_out(mpet_folder())
                            if spec_cap[-1] < 0.5*sp_cap():
                                print('simulation did not finish again, returning 0.25')
                                return 0.25
                        
                    elif prot == 'pitt':
                        time_tolerance = 3*60 # seconds
                        print('TIME END SIMULATION: ', tim[-1])
                        if tim[-1] < 1e3-time_tolerance:
                            print('simulation did not finish')
                            print('decreasing time tolerance and running again')
                            relTol, absTol = take_original_tolerance(config_files[0])
                            update_tolerance(ensambles_tot_new, config_files[0])
                            call_process(mpet_folder(), config_files[0])
                            reset_tolerance(ensambles_tot_new, config_files[0], relTol, absTol)
                            volt, curr, spec_cap, tim = take_data_sim_out(mpet_folder())
                            time_tolerance*=1.5
                            if tim[-1] < 1e3-time_tolerance:
                                print('simulation did not finish again, returning 0.25')
                                return 0.25
                    else:
                        print('protocol not recognized')

                if prot == 'cc' or prot == 'gitt':
                    ouput = volt
                elif prot == 'pitt':
                    ouput = curr
                dict_res[prot][temp][cond]["out"] = ouput
                dict_res[prot][temp][cond]["spec_cap"] = spec_cap
                dict_res[prot][temp][cond]["time"] = tim
                copy_sim_out_in_folder(mpet_folder(), store_folder)

    return dict_res


def obj(parameters):
    
    ensamble_tot, operating_conditions, config_files, optimization_method = sim_instructions()

    protocols = operating_conditions[0]
    c_rates_cc = operating_conditions[1]
    c_rates_gitt = operating_conditions[2]
    cont_volts = operating_conditions[3]
    temps = operating_conditions[4]


    dict_data = import_data(temps, protocols, c_rates_cc, c_rates_gitt, cont_volts)

    dict_res = run_mpet(parameters)
    if dict_res == 0.25:
        return 0.25
    
    weighted_mse = 0
    for prot in protocols:
        if prot == 'cc':
            applied_conditions = c_rates_cc
        elif prot == 'gitt':
            applied_conditions = c_rates_gitt
        elif prot == 'pitt':
            applied_conditions = cont_volts
        else:
            print('protocol not recognized')
        # header = "Temp(K)\tCondition(C-Rate/Volts)\tSpecCap(mAh/cm2)\tStepTime(s)\tOutput(V/C-rate)\n"

        for i, temp in enumerate(temps):
            for cond in applied_conditions[i]:
                print('calculating mse for protocol: ', prot, ' at temperature: ', temp, ' at condition: ', cond)
                out_sim = dict_res[prot][temp][cond]["out"]
                spcap_sim = dict_res[prot][temp][cond]["spec_cap"]
                spcap_sim -= c0()
                time_sim = dict_res[prot][temp][cond]["time"]

                out_exp = dict_data[prot][temp][cond]['Output(V/C-rate)']
                time_exp = dict_data[prot][temp][cond]['StepTime(s)'] 
                spcap_exp = dict_data[prot][temp][cond]['SpecCap(mAh/cm2)'] 
                    
                if prot == 'cc' or prot == 'gitt':
                    last_sim_cap = spcap_sim[-1]
                    if last_sim_cap < spcap_exp[-1]:
                        index_exp = np.where(spcap_exp>last_sim_cap)[0][0]
                        spcap_exp = spcap_exp[:index_exp]
                        time_exp = time_exp[:index_exp]
                        out_exp = out_exp[:index_exp]
                    out_interp = np.interp(spcap_exp, spcap_sim, out_sim)
                    # plt.plot(spcap_exp, out_exp, 'o',label='exp')
                    # plt.plot(spcap_exp, out_interp, 'o',label='sim')
                    # plt.plot(spcap_sim, out_sim, label='sim')
                    # plt.legend()
                    # plt.show()
                
                elif prot == 'pitt':
                    last_sim_time = time_sim[-1]
                    if last_sim_time < time_exp[-1]:
                        index_exp = np.where(time_exp>last_sim_time)[0][0]
                        spcap_exp = spcap_exp[:index_exp]
                        time_exp = time_exp[:index_exp]
                        out_exp = out_exp[:index_exp]
                    
                    out_interp = np.interp(time_exp, time_sim, out_sim)\
                    # remove initialization point

                    out_interp = out_interp[1:]
                    time_sim = time_sim[1:]
                    out_sim = out_sim[1:]
                    time_exp = time_exp[1:]
                    out_exp = out_exp[1:]

                    # plt.plot(time_exp, out_exp, 'o',label='exp')
                    # plt.plot(time_exp, out_interp, 'o',label='sim')
                    # plt.plot(time_sim, out_sim, label='sim')
                    # plt.legend()
                    # plt.show()
                else:
                    print('protocol not recognized')
                
                weighted_mse += np.average((out_exp - out_interp)**2) * weight_mse(prot, temp, cond)
            weighted_mse /= len(applied_conditions)
            
    weighted_mse /= len(protocols)*len(temps)
    root_weighted_mse = np.sqrt(weighted_mse)

    # weighted_mse = weighted_mse*10
    print('round done')
    with open(os.path.join(save_folder(),"log_book.txt"), 'a') as f:
        param_rescaled = copy.deepcopy(parameters)
        scaled_enasmble, scaling_factors = scale(ensamble_tot)
        params = rescale(param_rescaled, scaling_factors)
        for par in params:
            f.write(str(par) + "\t")
        f.write(str(root_weighted_mse) + "\n")

    return root_weighted_mse