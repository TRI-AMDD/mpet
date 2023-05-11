import os
import subprocess
import itertools
import configparser
# import numpy as np

def run_MPET(cwd, config):
    os.chdir(cwd)
    subprocess.call(["python", cwd + r"\bin\mpetrun.py", config])

def ensemble_definitions(parameters):
    ensamble = parameters
    keys = [vals[0] for vals in ensamble]
    val = [vals[1] for vals in ensamble]
    return keys, val

def run_params_mpet(config_file, material_file, system_properties, material_properties, output_folder):
    cwd = os.getcwd()
    os.chdir(cwd)
    if not os.path.exists(cwd + "\\" + output_folder):
        os.mkdir(cwd + "\\" + output_folder, )
    os.chdir(r".\configs")

    keys_system, val_system = ensemble_definitions(system_properties)
    cfg_sys = configparser.ConfigParser()
    cfg_sys.optionxform = str
    cfg_sys.read(config_file)
    combinations_system = list(itertools.product(*val_system))

    keys_mat, val_mat = ensemble_definitions(material_properties)
    cfg_mat = configparser.ConfigParser()
    cfg_mat.optionxform = str
    cfg_mat.read(material_file)
    combinations_material = list(itertools.product(*val_mat))
    for comb_mat in combinations_material:
        param_mat = dict(zip(keys_mat, comb_mat))
        new_mat = cfg_mat
        nicename_mat = []
        for key, val in param_mat.items():
            new_mat[key[0]][key[1]] = val
            nicename_mat.append(key[1] + "=" + val)
        with open(material_file, "w") as f:
            new_mat.write(f)
        for combin_sys in combinations_system:
            params_sys = dict(zip(keys_system, combin_sys))
            new_sys = cfg_sys
            nicename_sys = []
            for key, val in params_sys.items():
                new_sys[key[0]][key[1]] = val
                nicename_sys.append(key[1] + "=" + val)
            # cfg_dir = os.path.dirname(config_file)
            with open(config_file, "w") as f:
                new_sys.write(f)
    # nicename = update_config(config_file, parameters)
            os.chdir(cwd)
            run_MPET(cwd, (cwd + r'\\configs\\' + config_file))
        # change name of the last created folder in the folder "history"
            os.chdir(r".\history")
        # find the last created folder
            folders = os.listdir()
            folders.sort()
            last_folder = folders[0]
        # change the name of the last created folder
            new_folder_name = ""
            nicename = nicename_sys + nicename_mat
            for i in range(len(nicename)):
                new_folder_name = new_folder_name + nicename[i] + "_"
            os.rename(last_folder, os.path.join(cwd, output_folder, new_folder_name))
            os.chdir(cwd)
            os.chdir(r".\configs")


system_properties = [
        [("Conductivity", "sigma_s_c"), ["0.4","0.5"]],
        [("Sim Params","Crate"), ["1","2"]],
                ]

material_properties = [
    [("Reactions", 'k0'), ["0.1", "1"]],
    [("Reactions", 'E_A'), ["13000", "20000"]],
                ]

output_folder = "lfp_test"
config_file = 'params_system.cfg'
material_file = 'params_LFP.cfg'



run_params_mpet(config_file, material_file, system_properties, material_properties, output_folder)
