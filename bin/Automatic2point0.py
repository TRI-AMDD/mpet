import os
import subprocess
import itertools
import configparser
import time
# import numpy as np

def run_MPET(cwd, config):
    os.chdir(cwd)
    subprocess.call(["python", cwd + r"\bin\mpetrun.py", config])

def plot_data(output_folder):
    cwd = os.getcwd()
    plotter_folder = r"C:\Users\pierfrancescoo\Documents\Phase-field\Plot\plotter-TUD\Plotter"
    subprocess.call(["python", plotter_folder + r"\main.py", os.path.join(cwd, output_folder)])
    os.chdir(cwd)

def ensemble_definitions(parameters):
    ensamble = parameters
    keys = [vals[0] for vals in ensamble]
    val = [vals[1] for vals in ensamble]
    return keys, val


def run_params_mpet(config_file, material_file, battery_properties, material_properties, rate_properties, output_folder):
    cwd = os.getcwd()
    os.chdir(cwd)
    if not os.path.exists(cwd + "\\" + output_folder):
        os.mkdir(cwd + "\\" + output_folder, )
    os.chdir(r".\configs")

    keys_battery, val_battery = ensemble_definitions(battery_properties)
    cfg_bat = configparser.ConfigParser()
    cfg_bat.optionxform = str
    cfg_bat.read(config_file)
    combinations_battery = list(itertools.product(*val_battery))

    keys_mat, val_mat = ensemble_definitions(material_properties)
    cfg_mat = configparser.ConfigParser()
    cfg_mat.optionxform = str
    cfg_mat.read(material_file)
    combinations_material = list(itertools.product(*val_mat))

    keys_rate, val_rate = ensemble_definitions(rate_properties)
    cfg_rate = configparser.ConfigParser()
    cfg_rate.optionxform = str
    cfg_rate.read(config_file)
    combinations_rate = list(itertools.product(*val_rate))

    for comb_mat in combinations_material:
        param_mat = dict(zip(keys_mat, comb_mat))
        new_mat = cfg_mat
        nicename_mat = []
        nicename_bat = []
        for key, val in param_mat.items():
            new_mat[key[0]][key[1]] = val
            nicename_mat.append(key[1] + "=" + val)
        with open(material_file, "w") as f:
            new_mat.write(f)
        for combin_bat in combinations_battery:
            params_bat = dict(zip(keys_battery, combin_bat))
            new_bat = cfg_bat
            nicename_bat = []
            for key, val in params_bat.items():
                new_bat[key[0]][key[1]] = val
                nicename_bat.append(key[1] + "=" + val)
            with open(config_file, "w") as f:
                new_bat.write(f)
            # create folder of the simulation folders
            new_folder = ""
            nicename = []
            nicename = nicename_mat + nicename_bat
            
            for i in range(len(nicename)):
                new_folder = new_folder + nicename[i] + "-"
            if not os.path.exists(cwd + "\\" + output_folder + "\\" + new_folder):
                os.mkdir(cwd + "\\" + output_folder + "\\" + new_folder)

            for combin_rate in combinations_rate:
                params_rate = dict(zip(keys_rate, combin_rate))
                new_rate = cfg_rate
                for key, val in params_rate.items():
                    new_rate[key[0]][key[1]] = val
                    nicename_bat.append(key[1] + "=" + val)
                with open(config_file, "w") as f:
                    new_rate.write(f)
                # run MPET
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
                nicename = []
                nicename = nicename_mat + nicename_bat
                for i in range(len(nicename)):
                    new_folder_name = new_folder_name + nicename[i] + "-"
                if os.path.exists(os.path.join(cwd, output_folder + "\\" + new_folder, new_folder_name)):
                    os.rename(os.path.join(cwd, output_folder + "\\" + new_folder, new_folder_name),
                              os.path.join(cwd, output_folder + "\\" + new_folder, (new_folder_name + 'old')))
                os.rename(last_folder, os.path.join(cwd, output_folder + "\\" + new_folder, new_folder_name))
                os.chdir(cwd)
                os.chdir(r".\configs")

# polarization curve
battery_properties = [
        [("Conductivity","sigma_s_c"), ["0.1"]],
        [("Electrodes", "k0_foil"), ["14"]],
                ]

material_properties = [
        [("Reactions", 'k0'), ["50"]],
        [("Material", 'kappa'), ["20.0148e-10", "5.0148e-10"]],
                ]

rate_properties = [
        [("Sim Params","segments"), ["[(-0.2,150), (0,30), (-3,15)]",
                                     "[(-3,10), (0,30), (-3,15)]",
                                     "[(-5,6), (0,30), (-3,15)]"]],
                ]

# temperature_properties = [
#         [("Sim Params","T"), ["298"]],
#                 ]

output_folder = "LFP_memory1"
config_file = 'params_system_LFP_memory.cfg'
material_file = 'params_LFP_CHR_memory.cfg'



run_params_mpet(config_file, material_file, battery_properties,
                material_properties, rate_properties, output_folder)

# plot_data(output_folder)
