import os
import subprocess
import itertools
import configparser
import time


def run_MPET(cwd, config):
    os.chdir(cwd)
    subprocess.call(["python", os.path.join(cwd, "bin", "mpetrun.py"), config])

def ensemble_definitions(parameters):
    keys, vals = zip(*parameters)
    return keys, vals

def run_params_mpet(config_file, material_file, system_properties, material_properties, output_folder):
    cwd = os.getcwd()
    os.chdir(cwd)
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
    os.chdir("configs")

    keys_system, val_system = ensemble_definitions(system_properties)
    cfg_sys = configparser.ConfigParser()
    cfg_sys.optionxform = str
    cfg_sys.read(config_file)
    combinations_system = list(itertools.product(*val_system))
    num_sys = len(combinations_system)

    keys_mat, val_mat = ensemble_definitions(material_properties)
    cfg_mat = configparser.ConfigParser()
    cfg_mat.optionxform = str
    cfg_mat.read(material_file)
    combinations_material = list(itertools.product(*val_mat))
    num_mat = len(combinations_material)
    ind = 0

    for comb_mat in combinations_material:
        param_mat = dict(zip(keys_mat, comb_mat))
        new_mat = cfg_mat
        nicename_mat = []
        for key, val in param_mat.items():
            new_mat[key[0]][key[1]] = val
            nicename_mat.append(f"{key[1]}={val}")
        with open(material_file, "w") as f:
            new_mat.write(f)

        for combin_sys in combinations_system:
            params_sys = dict(zip(keys_system, combin_sys))
            new_sys = cfg_sys
            nicename_sys = []
            for key, val in params_sys.items():
                new_sys[key[0]][key[1]] = val
                nicename_sys.append(f"{key[1]}={val}")
            with open(config_file, "w") as f:
                new_sys.write(f)

            os.chdir(cwd)
            run_MPET(cwd, os.path.join(cwd, "configs", config_file))

            os.chdir("history")
            folders = os.listdir()
            folders.sort()
            last_folder = folders[0]
            new_folder_name = "-".join(nicename_mat + nicename_sys)
            new_folder_path = os.path.join(cwd, output_folder, new_folder_name)
            if os.path.exists(new_folder_path):
                os.rename(new_folder_path, os.path.join(cwd, output_folder, (new_folder_name + 'old')))
            os.rename(last_folder, new_folder_path)
            os.chdir(cwd)
            os.chdir("configs")
            ind += 1
            print(f"Simulation {ind} of {num_mat * num_sys} completed")



# system_properties = [
#         # [("Conductivity","G_mean_c"), ["0.5e-10"]],
#         # [("Conductivity","sigma_s_c"), ["1"]],
#         # [("Electrodes", "k0_foil"), ["20"]],
#         # [("Sim Params","Crate"), ["-5","-3"]],
#         # [("Particles","mean_c"), ["20e-9","50e-9","100e-9","300e-9"]],
#         # [("Particles","stddev_c"), ["35e-9"]],
#         # [("Geometry","BruggExp_c"), ["-1.85"]],
#         [("Sim Params","segments"), [
#                                     #  "[(-0.2,135), (0,60), (-3,10)]",
#                                      "[(-1,27), (0,60), (-3,10)]",
#                                     #  "[(-3,9), (0,60), (-3,10)]",
#                                     #  "[(-5,5.4), (0,60), (-3,10)]", 
#                                      ]],
#                 ]
 
# material_properties = [
#         # [("Reactions", 'k0'), ["30"]],
#         # [("Material", 'dgammadc'), ["0e-29"]],
#         # [("Material", 'kappa'), ["2.5148e-10"]],
#         # [("Material", 'D'), ["3.49e-15"]],
#         # [("Material", 'muRfunc'), ["LiFePO4_meta"]],
#         # [("Material", 'noise'), ["false"]],
#         # [("Material", 'noise_prefac'), ["1e-3"]],
#         # [("Material", 'size_dep_D'), ["true"]],
#         [("Material", 'B'), ["0.001e9"]],
#                 ]
 
# output_folder = "LFP_memory_best3C_1C"
# config_file = 'params_system_LFP_memory.cfg'
# material_file = 'params_LFP_CHR_memory.cfg'

system_properties = [
        [("Conductivity","G_mean_c"), ["4e-12","5e-12"]],
        # [("Conductivity","G_stddev_c", ["2e-13"]],
        # [("Conductivity","sigma_s_c"), ["2"]],
        # [("Electrodes", "k0_foil"), ["20"]],
        # [("Sim Params","Crate"), ["-3"]],
        # [("Particles","mean_c"), ["40e-9"]],
        # [("Particles","stddev_c"), ["50e-9"]],
        # [("Geometry","BruggExp_c"), ["-2"]],
        [("Sim Params","segments"), [
                                     "[(-5,5.4),(0,60),(-3,10)]", 
                                    #  "[(-3,9),(0,60),(-3,10)]",
                                     "[(-0.2,135),(0,60),(-3,10)]",
                                     ]],
                ]
 
material_properties = [
        [("Reactions", 'k0'), ["2.5","5"]], 
        # [("Material", 'dgammadc'), ["0e-29"]],
        # [("Material", 'kappa'), ["5e-10"]],
        [("Material", 'D_surf'), ["1e-18","1e-17"]],
        # [("Material", 'muRfunc'), ["LiFePO4_meta"]],
        # [("Material", 'noise'), ["false"]],
        # [("Material", 'noise_prefac'), ["1e-3"]],
        # [("Material", 'size_dep_D'), ["true"]], 
        [("Material", 'B'), ["0.005e9","0","0.02e9"]],
        [("Particles", 'thickness'), ["50e-9","100e-9"]],
                ]
 
output_folder = "LFP_ACR_noisek09_mem13"
config_file = 'params_system_LFP_Galuppini.cfg'
material_file = 'params_LFP_Galuppini.cfg'



run_params_mpet(config_file, material_file, system_properties, material_properties, output_folder)

# plot_data(output_folder)
