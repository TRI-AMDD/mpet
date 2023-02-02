import os
import subprocess


def run_MPET(cwd, config):
    os.chdir(cwd)
    subprocess.call(["python", cwd + r"\bin\mpetrun.py", config])


cwd = os.getcwd()
os.chdir(cwd)


# os.chdir(r".\configs\configs_memory")

conf = open('ensemble_parallel_configs.txt',"r")
lines = conf.readlines()
os.chdir(cwd)
# print(cwd)
for line in lines:
    cwd = os.getcwd()
    run_MPET(cwd, (cwd + r'\\' + line[:-1]))



            
