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

def run_params_mpet(config_file, parameters, output_folder):
    cwd = os.getcwd()
    os.chdir(cwd)
    if not os.path.exists(cwd + "\\" + output_folder):
        os.mkdir(cwd + "\\" + output_folder, )
    os.chdir(r".\configs")
    with open(config_file, "r") as ff:
        keys, val = ensemble_definitions(parameters)
        cfg = configparser.ConfigParser()
        cfg.optionxform = str
        cfg.read(config_file)
        combinations = list(itertools.product(*val))
        for combination in combinations:
            params = dict(zip(keys,combination))
            new_cfg = cfg
            nicename = []
            for key, val in params.items():
                new_cfg[key[0]][key[1]] = val
                nicename.append(key[1] + "=" + val)
            cfg_dir = os.path.dirname(config_file)
            with open(config_file, "w") as f:
                new_cfg.write(f)
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
            new_folder_name = str(nicename[0] + "_" + nicename[1])
            os.rename(last_folder, os.path.join(cwd, output_folder, new_folder_name))
            os.chdir(cwd)
            os.chdir(r".\configs")


parameters = [
        [("Conductivity", "sigma_s_c"), ["0.1","0.2","0.3", "0.4","0.5","0.6","0.7","0.8","0.9","1"]],
        [("Sim Params","Crate"), ["-0.62"]],
                ]

output_folder = "1C_sigmas"
config_file = 'params_system_NMC_Mark.cfg'



run_params_mpet(config_file ,parameters, output_folder)
