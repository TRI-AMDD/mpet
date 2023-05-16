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
        nicename_sys = []
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
            with open(config_file, "w") as f:
                new_sys.write(f)
                    
            os.chdir(cwd)
            # suspend for 20 seconds
            # time.sleep(0.5)

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
            nicename = nicename_mat + nicename_sys
            for i in range(len(nicename)):
                new_folder_name = new_folder_name + nicename[i] + "-"
            if os.path.exists(os.path.join(cwd, output_folder, new_folder_name)):
                os.rename(os.path.join(cwd, output_folder, new_folder_name), os.path.join(cwd, output_folder, (new_folder_name + 'old')))
            os.rename(last_folder, os.path.join(cwd, output_folder, new_folder_name))
            os.chdir(cwd)
            os.chdir(r".\configs")

# NMC study
# system_properties = [
#         [("Conductivity", "sigma_s_c"), ["0.06","0.05", "0.04"]],
#         [("Sim Params","Crate"), ["-0.206"]],
#         [("Electrodes", "Rfilm_foil"), ["0"]],
#         [("Electrodes", "k0_foil"), ["4","3","2"]],
#                 ]

# material_properties = [
#     [("Reactions", 'k0'), ["25","20","15"]],
#                 ]

# output_folder = "NMC_thick_conv_8032522"
# config_file = 'params_system_NMC_Mark.cfg'
# material_file = 'params_NMC_Chen2020.cfg'

# LFP study
# diffusion more or less 1e-15
# k0 between 40 and 70
# sigma_c at least 0.05, small influence
# system_properties = [
#         [("Sim Params","segments"), ["[(-1,60)]","[(-3,20)]"]],
#         [("Conductivity","sigma_s_c"), ["0.1","0.075", "0.05"]],
#                 ]

# material_properties = [
#     [("Reactions", 'k0'), ["40", "50", "70"]],
#     [("Material", 'D'), ["2e-15","0.5e-15", "1e-15", "3e-15"]],
#                 ]

# output_folder = "LFP_test_more_precise"
# config_file = 'params_system_LFP_test.cfg'
# material_file = 'params_LFP_CHR.cfg'
# k0_foil is  14
# k0 is 50
# sigma_s_c = 0.1
# polarization curve


system_properties = [
        # [("Conductivity","sigma_s_c"), ["0.1"]],
        [("Electrodes", "k0_foil"), ["15"]],
        [("Sim Params","segments"), ["[(-1,60)]"]],
                ]

material_properties = [
        # [("Reactions", 'k0'), ["25"]],
        # [("Material", 'kappa'), ["20.0148e-10", "5.0148e-10"]],
        [("Material", 'D'), ["5e-14"]],
                ]

output_folder = "LFP_memory10"
config_file = 'params_system_LFP.cfg'
material_file = 'params_LFP_CHR.cfg'



run_params_mpet(config_file, material_file, system_properties, material_properties, output_folder)

# plot_data(output_folder)
