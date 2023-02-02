import os
import subprocess


def run_MPET(cwd, config):
    os.chdir(cwd)
    subprocess.call(["python", cwd + "\mpet-dev-delft\\bin\mpetrun.py", config])


currentDir = os.getcwd()


cwd = os.getcwd()
os.chdir(cwd)
os.chdir(cwd+"\mpet-dev-delft\configs\configs_memory")
conf = open('configs.txt',"r")
lines = conf.readlines()
os.chdir(cwd)
# print(cwd)
for line in lines:
    print(line[:-1])
    run_MPET(cwd, (cwd + '\mpet-dev-delft\configs\configs_400p_stdCHR\\' + line[:-1]))



            
