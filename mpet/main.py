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

import mpet
import mpet.data_reporting as data_reporting
import mpet.maccor_reader as maccor_reader
import mpet.io_utils as IO
import mpet.sim as sim
import mpet.utils as utils


def run_simulation(ndD_s, ndD_e, tScale, outdir):
    # Create Log, Solver, DataReporter and Simulation object
    log = dae.daePythonStdOutLog()
    daesolver = dae.daeIDAS()
    simulation = sim.SimMPET(ndD_s, ndD_e, tScale)
    datareporter = data_reporting.setup_data_reporters(simulation, ndD_s, outdir)

    # Use SuperLU direct sparse LA solver
    lasolver = pySuperLU.daeCreateSuperLUSolver()
#    lasolver = pyTrilinos.daeCreateTrilinosSolver("Amesos_Umfpack", "")
    daesolver.SetLASolver(lasolver)

    # Enable reporting of all variables
    simulation.m.SetReportingOn(True)

    #Turn off reporting of some variables
    simulation.m.endCondition.ReportingOn=False
    
    #Turn off reporting of particle ports
    for trode in simulation.m.trodes:
        for particle in simulation.m.particles[trode]:
            pModel = particle[0]
            for port in pModel.Ports:
                for var in port.Variables:
                    var.ReportingOn = False
    
    #Turn off reporting of cell ports
    for port in simulation.m.Ports:
        for var in port.Variables:
            var.ReportingOn = False

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

    out = True #returns true if run was successful

    # Run
    try:
        simulation.Run()
    except Exception as e:
        print(str(e))
        simulation.ReportData(simulation.CurrentTime)
        #raise, commented out so that the system prints output even when
        #it errors out
        out = False
    except KeyboardInterrupt:
        print("\nphi_applied at ctrl-C:",
              simulation.m.phi_applied.GetValue(), "\n")
        simulation.ReportData(simulation.CurrentTime)
        out = False
    simulation.Finalize()

    return out

def main(paramfile, keepArchive=True):
    timeStart = time.time()
    # Get the parameters dictionary (and the config instance) from the
    # parameter file
    P_s, P_e = IO.get_configs(paramfile)
    dD_s, ndD_s, dD_e, ndD_e = IO.get_dicts_from_configs(P_s, P_e, paramfile)

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

    fo=open(os.path.join(outdir, 'run_info.txt'), 'w')
    
    #Print mpet version
    print("mpet version:", file=fo)
    print(mpet.__version__+"\n", file=fo)

    #Print git commit info if it exists
    if commit_hash != "":
        # Store commit info to file, as well as how to patch if
        # there's a diff
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
    
    fo.close()
    
    # External functions are not supported by the Compute Stack approach.
    # Activate the Evaluation Tree approach if noise, logPad, CCsegments,
    # or CVsegments are used
    cfg = dae.daeGetConfig()
    noise=ndD_e['c']['noise']
    logPad=ndD_e['c']['logPad']
    segments = ndD_s["profileType"] in ["CCsegments","CVsegments","CCCVCPcycle"]
    if (noise or logPad or (segments and ndD_s["tramp"]>0)) and cfg.has_key('daetools.core.equations.evaluationMode'):
        cfg.SetString('daetools.core.equations.evaluationMode', 'evaluationTree_OpenMP')
    
    cfg.SetInteger('daetools.IDAS.numberOfSTNRebuildsDuringInitialization', 10000)
    with open(os.path.join(outdir, "daetools_config_options.txt"), 'w') as fo:
        print(cfg, file=fo)

    P_s.set('Electrodes', 'cathode', 'input_params_c.cfg')
    P_s.set('Electrodes', 'anode', 'input_params_a.cfg')

    #Disable printStats
    cfg.SetString('daetools.activity.printStats','false')

    #save our input files and runs simulations
    if ndD_s["profileType"][-4:] != ".000":
        #if it is not a maccor cycling procedure file
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
        #carry out simulation
        run_simulation(ndD_s, ndD_e, dD_s["td"], outdir)
    else:
        #if it is a maccor cycling procedure file
        cycling_dicts = maccor_reader.get_cycling_dict(ndD_s, dD_s)
        for i in range(len(cycling_dicts)):
            if i == 0:
                #only saves data for electrodes once since it doesnt change
                for trode in ndD_s["trodes"]:
                    paramFileName = "input_params_{t}.cfg".format(t=trode)
                    paramFile = os.path.join(outdir, paramFileName)
                    IO.write_config_file(P_e[trode], filename=paramFile)
                    dictFile = os.path.join(outdir, "input_dict_{t}".format(t=trode))
                    IO.write_dicts(dD_e[trode], ndD_e[trode], filenamebase=dictFile)       
            else: #if its not the first simulation, then set prevDir
                P_s.set('Sim Params', 'prevDir', outdir)
            num_steps = len(cycling_dicts["step_" + str(i)].get("segments"))
            P_s.set('Sim Params', 'period', str([1e80]*num_steps))
            #sets the segments and total cycle numbers for each set we are going to run
            P_s.set('Sim Params', 'segments', str(cycling_dicts["step_" + str(i)].get("segments")))
            P_s.set('Sim Params', 'totalCycle', cycling_dicts["step_" + str(i)].get("totalCycle"))
            P_s.set('Sim Params', 'profileType', 'CCCVCPcycle')
            #fills in period for waveform
            dD_s, ndD_s, dD_e, ndD_e = IO.get_dicts_from_configs(P_s, P_e)
            #reset everything in dictionaries too
            paramFileName = "input_params_system_{j}.cfg".format(j=str(i))
            paramFile = os.path.join(outdir, paramFileName)
            IO.write_config_file(P_s, filename=paramFile)
            dictFile = os.path.join(outdir, "input_dict_system_" + str(i))
            IO.write_dicts(dD_s, ndD_s, filenamebase=dictFile)
            #carries out simulation
            out = run_simulation(ndD_s, ndD_e, dD_s["td"], outdir)
            if out == False:
                break

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
