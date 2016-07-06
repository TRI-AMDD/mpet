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

import mpet.data_reporting as data_reporting
import mpet.io_utils as IO
import mpet.sim as sim


def consoleRun(ndD_s, ndD_e, tScale, outdir):
    # Create Log, Solver, DataReporter and Simulation object
    log = dae.daePythonStdOutLog()
    daesolver = dae.daeIDAS()
    simulation = sim.SimMPET(ndD_s, ndD_e, tScale)
    datareporter = data_reporting.setupDataReporters(simulation, outdir)

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
    simulation.ReportingInterval = ndD_s["tend"] / ndD_s['tsteps']

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


def main(paramfile="params_default.cfg", keepArchive=True):
    timeStart = time.time()
    # Get the parameters dictionary (and the config instance) from the
    # parameter file
    P_s, P_e = IO.getConfigs(paramfile)
    dD_s, ndD_s, dD_e, ndD_e = IO.getDictsFromConfigs(P_s, P_e)

    # Directories we'll store output in.
    outdir_name = time.strftime("%Y%m%d_%H%M%S", time.localtime())
    outdir_path = os.path.join(os.getcwd(), "history")
    outdir = os.path.join(outdir_path, outdir_name)
    # Make sure there's a place to store the output
    try:
        os.makedirs(outdir)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
        # Clean out the directory
        for file_obj in os.listdir(outdir):
            file_path = os.path.join(outdir, file_obj)
            if os.path.isfile(file_path):
                os.unlink(file_path)
            else:
                shutil.rmtree(file_path)
    paramFileName = "input_params_system.cfg"
    paramFile = os.path.join(outdir, paramFileName)
    IO.writeConfigFile(P_s, filename=paramFile)
    dictFile = os.path.join(outdir, "input_dict_system")
    IO.writeDicts(dD_s, ndD_s, filenamebase=dictFile)
    for trode in ndD_s["trodes"]:
        paramFileName = "input_params_{t}.cfg".format(t=trode)
        paramFile = os.path.join(outdir, paramFileName)
        IO.writeConfigFile(P_e[trode], filename=paramFile)
        dictFile = os.path.join(outdir, "input_dict_{t}".format(t=trode))
        IO.writeDicts(dD_e[trode], ndD_e[trode], filenamebase=dictFile)

    # Store info about this script
    # mpet.py script directory
    localDir = os.path.dirname(os.path.abspath(__file__))
    out1 = ""
    try:
        # Git option, if it works -- commit info and current diff
        out1 = subp.check_output(['git', '-C', localDir, 'rev-parse', '--short', 'HEAD'],
                                 stderr=subp.STDOUT, universal_newlines=True)
        out2 = subp.check_output(['git', '-C', localDir, 'diff'],
                                 stderr=subp.STDOUT, universal_newlines=True)
        out3 = subp.check_output(
            ['git', '-C', localDir, 'rev-parse', '--abbrev-ref', 'HEAD'],
            stderr=subp.STDOUT, universal_newlines=True)
    except subp.CalledProcessError:
        pass
    if out1 != "":
        # Store commit info to file, as well as how to patch if
        # there's a diff
        with open(os.path.join(outdir, 'run_info.txt'), 'w') as fo:
            print("branch name:", file=fo)
            print(out3, file=fo)
            print("commit hash:", file=fo)
            print(out1, file=fo)
            print("to run:", file=fo)
            print("$ git checkout [commit hash]", file=fo)
            print("$ patch -p1 < commit.diff:", file=fo)
            print("$ python[3] mpet.py input_params.cfg", file=fo)
        with open(os.path.join(outdir, 'commit.diff'), 'w') as fo:
            print(out2, file=fo)
    else:
        # At least keep a copy of the python files in this directory
        # with the output
        snapshotDir = os.path.join(outdir, "simSnapshot")
        os.makedirs(snapshotDir)
        pyFiles = glob.glob(os.path.join(localDir, "*.py"))
        for pyFile in pyFiles:
            shutil.copy(pyFile, snapshotDir)
    if sys.platform in ["linux", "linux2"]:
        cfgLoc = os.path.join("/", "etc", "daetools", "daetools.cfg")
    elif sys.platform in ["win32"]:
        cfgLoc = os.path.join("/", "daetools", "daetools.cfg")
    elif sys.platform in ["cygwin"]:
        cfgLoc = os.path.join("/", "cygdrive", "c", "daetools", "daetools.cfg")
    try:
        shutil.copy(cfgLoc, outdir)
    except:
        if sys.platform in ["linux", "linux2"]:
            try:
                cfgdir = os.path.join(os.environ["HOME"], ".daetools")
                shutil.copy(os.path.join(cfgdir, "daetools.cfg"), outdir)
            except:
                pass
        else:
            pass

    # Carry out the simulation
    consoleRun(ndD_s, ndD_e, dD_s["td"], outdir)

    # Final output for user
    if paramfile == "params_default.cfg":
        print("\n\n*** WARNING: Used default file, ""{fname}"" ***".format(
            fname=default_file))
        print("Pass other parameter file as an argument to this script\n")
    else:
        print("\n\nUsed parameter file ""{fname}""\n\n".format(
            fname=paramfile))
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
        shutil.copy(os.path.join(outdir, fname), tmpDir)

    if not keepArchive:
        shutil.rmtree(outdir)
