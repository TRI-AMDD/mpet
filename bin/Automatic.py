import os
import subprocess
import itertools
# import numpy as np


def run_MPET(cwd, config):
    os.chdir(cwd)
    subprocess.call(["python", cwd + r"\bin\mpetrun.py", config])


def update_config(config_file, parameters):
    with open(config_file, "r") as f:
        lines = f.readlines()

    for param_name, param_value in parameters.items():
        for i, line in enumerate(lines):
            if line.startswith(f"{param_name} = "):
                lines[i] = f"{param_name} = {param_value}\n"
                break

    with open(config_file, "w") as f:
        f.writelines(lines)


cwd = os.getcwd()
os.chdir(cwd)

parameters = {
    "Crate": [-0.62, -0.206],
    "sigma_s_c": [0.1, 0.2, 0.3, 0.4, 0.5],
}

# parameters = [
#         [("Conductivity", "sigma_s_c"), ["0.1","0.2","0.3"]],
#         [("Geometry","BruggExp_c"), ["-0.05","-0.1","-0.15"]],
#                 ]
keys = [vals[0] for vals in parameters]
val = [vals[1] for vals in parameters]

combinations = list(itertools.product(*val))
print(combinations)

output_folder = "1C_sigmas"
config_file = 'params_system_NMC_Mark.cfg'
# mkdir even if it already exists
if not os.path.exists(cwd + "\\" + output_folder):
    os.mkdir(cwd + "\\" + output_folder, )

os.chdir(r".\configs")

for param_name, param_values in parameters.items():
    for param_value in param_values:
        update_config(config_file, {param_name: param_value})
        os.chdir(cwd)
        run_MPET(cwd, (cwd + r'\\configs\\' + config_file))
        # change name of the last created folder in the folder "history"
        os.chdir(r".\history")
        # find the last created folder
        folders = os.listdir()
        folders.sort()
        last_folder = folders[0]
        # change the name of the last created folder
        new_folder_name = output_folder + f"_{param_name}={param_value}"
        os.rename(last_folder, os.path.join(cwd, output_folder, new_folder_name))
        os.chdir(cwd)
        os.chdir(r".\configs")
