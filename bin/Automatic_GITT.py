import os
import subprocess
import itertools
import configparser


def run_MPET(cwd, config):
    os.chdir(cwd)
    subprocess.call(["python", os.path.join(cwd, "bin", "mpetrun.py"), config])


def ensemble_definitions(parameters):
    keys, vals = zip(*parameters)
    return keys, vals


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
                if key[1] == "segments":
                    nicename_sys.append(f"{key[1]}={val[9:13]}")
                else:
                    nicename_sys.append(f"{key[1]}={val}")
                new_sys["Sim Params"]["profileType"] = "CCsegments"
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

# Even in generous situations (high numb of cont) G = 1e-14 is too low, G = 1e-13 seems okish

n_pulses = 20
rest_time = 30  # min
Crates = [1]
normaliz = 160/170


list_segments = []
for Crate in Crates:
    time_per_pulse = 60 / (Crate * n_pulses)
    segments = f"[(0,{rest_time}),"
    for i in range(n_pulses+2):
        string = f"({Crate*normaliz},{time_per_pulse}),(0,{rest_time}),"
        segments += string  # Correct way to concatenate the string

    # Adding the final closing bracket
    segments += "]"
    list_segments.append(segments)


system_properties = [
    [("Sim Params","segments"), list_segments],
    # [("Conductivity","avg_num_cont_c"), ["2"]],
    # [("Conductivity","std_num_cont_c"), ["4"]],
    # [("Conductivity","avg_num_cont_c"), ["2"]],
    # [("Conductivity","std_num_cont_c"), ["2"]],
    # [("Conductivity","perc_grid_c"), ["0.3"]],
    # [("Conductivity","sigma_s_c"), ["1"]],
    # [("Conductivity","sig_carb_c"), ["1e-4"]],
    # [("Conductivity","sig_bulk_c"), ["5e-8"]],
    [("Conductivity","avg_num_cont_c"), ["2"]],
    [("Conductivity","std_num_cont_c"), ["2"]],
    [("Conductivity","c_dep_exp_c"), ["1"]],
    # [("Conductivity","perc_grid_c"), ["0.4","0.8"]],
    # [("Conductivity","sigma_s_c"), ["1"]],
    [("Conductivity","sig_carb_c"), ["2e-4"]],
    [("Conductivity","sig_bulk_c"), ["1e-6"]],
    [('Geometry',"L_c"), ["40e-6"]],
    # [("Sim Params","seed"), ["0"]],
    # [("Conductivity","sigma_s_c"), ["0.5"]],
    # [("Sim Params","Npart_c"), ["20"]],
    # [("Sim Params","Nvol_c"), ["2"]],
    # [("Particles","mean_c"), ["50e-9","100e-9"]],
    # [("Particles","stddev_c"), ["15e-9"]],
    # [("Sim Params","Crate"), ["1.5"]],
    # [("Sim Params","T"), ["268","283","298"]],
    ]

material_properties = [
    [("Reactions", 'k0'), ["15"]],
    # [("Material", 'D'), ["1e-17","5e-18"]],
    # [("Reactions", 'Rfilm'), ["0"]],
    ]


output_folder = "LFP_GITT/conc_dep_cond"
config_file = 'params_system_Temp.cfg'
material_file = 'params_LFP.cfg'


run_params_mpet(config_file, material_file, system_properties, material_properties, output_folder)
