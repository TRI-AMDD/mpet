"""The main module that organizes the simulation and manages data IO."""
import errno
import glob
import os
import shutil
import subprocess as subp
import sys
import time

import daetools.pyDAE as dae
from daetools.solvers.superlu import pySuperLU
import numpy as np

import mpet.data_reporting as data_reporting
import mpet.io_utils as IO
import mpet.sim as sim
import mpet.utils as utils


def run_simulation(ndD_s, ndD_e, tScale, outdir):
    # Create Log, Solver, DataReporter and Simulation object
    log = dae.daePythonStdOutLog()
    daesolver = dae.daeIDAS()
    simulation = sim.SimMPET(ndD_s, ndD_e, tScale)
    datareporter = data_reporting.setup_data_reporters(simulation, outdir)

    # Use SuperLU direct sparse LA solver
    lasolver = pySuperLU.daeCreateSuperLUSolver()
#    lasolver = pyTrilinos.daeCreateTrilinosSolver("Amesos_Umfpack", "")
    daesolver.SetLASolver(lasolver)

    # Enable reporting of all variables
    simulation.m.SetReportingOn(True)

    # Set relative tolerances
    daesolver.RelativeTolerance = ndD_s["relTol"]

    # Set the time horizon and the reporting interval
    simulation.TimeHorizon = ndD_s["tend"]
    # The list of reporting times excludes the first index (zero, which is implied)
    simulation.ReportingTimes = list(np.linspace(0, ndD_s["tend"], ndD_s["tsteps"] + 1))[1:]
    # Example logspacing for output times:
    # simulation.ReportingTimes = list(
    #     np.logspace(-4, np.log10(simulation.TimeHorizon), ndD_s['tsteps']))

    # Connect data reporter
    simName = simulation.m.Name + time.strftime(
        " [%d.%m.%Y %H:%M:%S]", time.localtime())
    if not datareporter.Connect("", simName):
        sys.exit()

    # Initialize the simulation
    simulation.Initialize(daesolver, datareporter, log)

    # Solve at time=0 (initialization)
    simulation.SolveInitial()

    # Run
    try:
        simulation.Run()
    except Exception as e:
        print(str(e))
        simulation.ReportData(simulation.CurrentTime)
    except KeyboardInterrupt:
        print("\nphi_applied at ctrl-C:",
              simulation.m.phi_applied.GetValue(), "\n")
        simulation.ReportData(simulation.CurrentTime)
    simulation.Finalize()


def main(paramfile, keepArchive=True):
    timeStart = time.time()
    # Get the parameters dictionary (and the config instance) from the
    # parameter file
    P_s, P_e = IO.get_configs(paramfile)
    dD_s, ndD_s, dD_e, ndD_e = IO.get_dicts_from_configs(P_s, P_e)

    # Directories we'll store output in.
    outdir_name = time.strftime("%Y%m%d_%H%M%S", time.localtime())
    outdir_path = os.path.join(os.getcwd(), "history")
    outdir = os.path.join(outdir_path, outdir_name)
    # Make sure there's a place to store the output
    try:
        os.makedirs(outdir)
    except OSError as exception:
        if exception.errno == errno.EEXIST:
            print("The output directory, {dirname}, exists. Aborting.".format(dirname=outdir))
            sys.exit()
        else:
            raise
    paramFileName = "input_params_system.cfg"
    paramFile = os.path.join(outdir, paramFileName)
    IO.write_config_file(P_s, filename=paramFile)
    dictFile = os.path.join(outdir, "input_dict_system")
    IO.write_dicts(dD_s, ndD_s, filenamebase=dictFile)
    for trode in ndD_s["trodes"]:
        paramFileName = "input_params_{t}.cfg".format(t=trode)
        paramFile = os.path.join(outdir, paramFileName)
        IO.write_config_file(P_e[trode], filename=paramFile)
        dictFile = os.path.join(outdir, "input_dict_{t}".format(t=trode))
        IO.write_dicts(dD_e[trode], ndD_e[trode], filenamebase=dictFile)

    # Store info about this script
    # mpet.py script directory
    localDir = os.path.dirname(os.path.abspath(__file__))
    commit_hash = ""
    try:
        # Git option, if it works -- commit info and current diff
        branch_name, commit_hash, commit_diff = utils.get_git_info(localDir, shell=False)
    except FileNotFoundError:
        try:
            branch_name, commit_hash, commit_diff = utils.get_git_info(localDir, shell=True)
        except subp.CalledProcessError:
            pass
    except subp.CalledProcessError:
        pass
    if commit_hash != "":
        # Store commit info to file, as well as how to patch if
        # there's a diff
        with open(os.path.join(outdir, 'run_info.txt'), 'w') as fo:
            print("branch name:", file=fo)
            print(branch_name, file=fo)
            print("commit hash:", file=fo)
            print(commit_hash, file=fo)
            print("to run, from the root repo directory, copy relevant files there,", file=fo)
            print("edit input_params_system.cfg to point to correct material", file=fo)
            print("params files, and:", file=fo)
            print("$ git checkout [commit hash]", file=fo)
            print("$ patch -p1 < commit.diff:", file=fo)
            print("$ python[3] mpetrun.py input_params_system.cfg", file=fo)
        with open(os.path.join(outdir, 'commit.diff'), 'w') as fo:
            print(commit_diff, file=fo)
    else:
        # At least keep a copy of the python files in this directory
        # with the output
        snapshotDir = os.path.join(outdir, "simSnapshot")
        os.makedirs(snapshotDir)
        pyFiles = glob.glob(os.path.join(localDir, "*.py"))
        for pyFile in pyFiles:
            shutil.copy(pyFile, snapshotDir)
    cfg = dae.daeGetConfig()
    with open(os.path.join(outdir, "daetools_config_options.txt"), 'w') as fo:
        print(cfg, file=fo)

    # Carry out the simulation
    run_simulation(ndD_s, ndD_e, dD_s["td"], outdir)

    # Final output for user
    print("\n\nUsed parameter file ""{fname}""\n\n".format(fname=paramfile))
    timeEnd = time.time()
    tTot = timeEnd - timeStart
    print("Total time:", tTot, "s")
    try:
        with open(os.path.join(outdir, 'run_info.txt'), 'a') as fo:
            print("\nTotal run time:", tTot, "s", file=fo)
    except Exception:
        pass

    # Copy simulation output to current directory
    tmpDir_name = "sim_output"
    tmpDir = os.path.join(os.getcwd(), tmpDir_name)
    try:
        os.makedirs(tmpDir)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
    for fname in os.listdir(outdir):
        tocopy = os.path.join(outdir, fname)
        if os.path.isfile(tocopy):
            shutil.copy(tocopy, tmpDir)

    if not keepArchive:
        shutil.rmtree(outdir)
