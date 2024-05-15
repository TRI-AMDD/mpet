import os
import subprocess
import itertools
import configparser
import numpy as np


def run_MPET(cwd, config):
    os.chdir(cwd)
    subprocess.call(["python", os.path.join(cwd, "bin", "mpetrun.py"), config])


def ensemble_definitions(parameters):
    keys, vals = zip(*parameters)
    return keys, vals

k_B = 1.38064852e-23
T = 298.15
e = 1.60217662e-19

def run_params_mpet(config_file, material_file,
                    system_properties, material_properties, output_folder):
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
                if key[1] == "prevDir":
                    continue
                new_sys["Sim Params"]["profileType"] = "CVsegments"
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
                os.rename(new_folder_path, os.path.join(cwd, output_folder,
                                                        (new_folder_name + 'old')))
            os.rename(last_folder, new_folder_path)
            os.chdir(cwd)
            os.chdir("configs")
            ind += 1
            print(f"Simulation {ind} of {num_mat * num_sys} completed")

# prev_dir = r"C:\Users\pierfrancescoo\Documents\Phase-field\mpet-LFMP\mpet\LFMP_dyn\pulses_y04\50percMn\base_50pMn"
# prev_dir = r"C:\Users\pierfrancescoo\Documents\Phase-field\mpet-LFMP\mpet\LFP_CV\Iarchuk_1\base"
# ocv = 3.9998
ocv = 3.422
etas = [2,3,4]

holds = 30000 # sec
holds = holds/60 # min
etas = k_B*T/e*np.array(etas)
segments = []
for i in range(len(etas)):
    Vp = str(ocv + etas[i])
    Vm = str(ocv - etas[i])
    holds = str(holds)
    # stringp = f"[(3.20,120),(3.4,30),({Vp},{holds})]"
    # stringm = f"[(3.47,20),({Vm},{holds})]"
    # stringp = f"[({Vp},{holds})]"
    stringm = f"[({Vm},{holds})]"
    # segments.append(str(stringp))
    segments.append(str(stringm))


system_properties = [
    
    # [("Sim Params","seed"), ["0"]],
    # [("Conductivity","perc_grid_c"), ["0.45"]],
    # [("Conductivity","avg_num_cont_c"), ["2"]],
    # [("Conductivity","sig_bulk_c"), ["3e-8"]],
    # [("Conductivity","sigma_s_c"), ["0.05"]],
    # [("Conductivity","avg_num_cont_c"), ["2"]],
    # [("Conductivity","std_num_cont_c"), ["2"]],
    # [("Conductivity","perc_grid_c"), ["0.3"]],
    [("Conductivity","c_dep_exp_c"), ["0"]],
    [("Conductivity","sig_carb_c"), ["1e-4"]],
    [("Conductivity","sig_bulk_c"), ["5e-8"]],
    [('Geometry',"L_c"), ["40e-6"]],
    # [("Sim Params","Npart_c"), ["5"]],
    # [("Sim Params","Nvol_c"), ["10"]],
    # [("Particles","mean_c"), ["100e-9"]],
    # [("Particles","stddev_c"), ["25e-9"]],
    # [("Sim Params","prevDir"), [prev_dir]],
    # [("Electrolyte","c0"), ["1000"]],
    [("Sim Params","segments"), segments],
    ]

material_properties = [
    [("Reactions", 'k0'), ["15"]],
    # [("Particles", 'std_thickness'), ["20e-9"]],
    # [("Material", 'delta_gamma_vert'), ["0e-30"]],
    # [("Reactions", 'Rfilm'), ["0"]],
    # [("Material", 'B'), ["0.1916e9"]],
    # [("Material", 'D'), ["1e-17"]],
    # [("Material", 'kappa'), ["5e-10"]],
    # [("Reactions", 'lambda'), ["5.54e-20"]],
    # [("Reactions", 'lambda'), ["3.4113e-20"]],
    # [("Material", 'cwet'), ["0.98"]],
    ]


output_folder = "LFP_PITT\\con_dep_cond"
config_file = 'params_system_Temp.cfg'
material_file = 'params_LFP.cfg'

run_params_mpet(config_file, material_file, system_properties, material_properties, output_folder)
