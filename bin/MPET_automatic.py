import os
import subprocess
import numpy as np


def run_MPET(cwd, config):
    os.chdir(cwd)
    subprocess.call(["python", cwd + r"\bin\mpetrun.py", config])


cwd = os.getcwd()
os.chdir(cwd)


os.chdir(r".\configs\mark_diff_Crate")
config_list = os.listdir()
# print(config_list)
config_I_want = np.array([])
for config in config_list:
    if config[:22] == 'params_system_NMC_Mark':
        config_I_want = np.append(config_I_want,config)
print(config_I_want)
# conf = open('ensemble_parallel_configs.txt',"r")
# lines = conf.readlines()
os.chdir(cwd)
# # print(cwd)
for line in config_I_want:
    cwd = os.getcwd()
    run_MPET(cwd, (cwd + r'\\' + r'configs\\mark_diff_Crate\\' + line))



            
