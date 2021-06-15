import errno
import os
import os.path as osp
import shutil
import time

import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio

import mpet.main
from mpet.config.configuration import Config
import tests.test_defs as defs


def run_test_sims(runInfo, dirDict, pflag=True):
    for testStr in runInfo:
        testDir = osp.join(dirDict["out"], testStr)
        os.makedirs(testDir)

        # Copy config files from ref_outputs
        refDir = osp.join(dirDict["suite"],"ref_outputs",testStr)
        _, _, filenames = next(os.walk(osp.join(refDir)))
        for f in filenames:
            if ".cfg" in f:
                shutil.copyfile(osp.join(refDir,f),osp.join(testDir,f))

        # Run the simulation
        configfile = osp.join(testDir,'params_system.cfg')
        mpet.main.main(configfile, keepArchive=False)
        shutil.move(dirDict["simOut"], testDir)

    # Remove the history directory that mpet creates.
    try:
        os.rmdir(osp.join(dirDict["suite"], "history"))
    except OSError as exception:
        if exception.errno != errno.ENOENT:
            raise


def run_test_sims_analyt(runInfo, dirDict, pflag=True):
    for testStr in sorted(runInfo.keys()):
        testDir = osp.join(dirDict["out"], testStr)
        os.makedirs(testDir)
        runInfo[testStr][0](testDir, dirDict)
    # Remove the history directory that mpet creates.
    try:
        os.rmdir(osp.join(dirDict["suite"], "history"))
    except OSError as exception:
        if exception.errno != errno.ENOENT:
            raise


def get_sim_time(simDir):
    with open(osp.join(simDir, "run_info.txt")) as fi:
        simTime = float(fi.readlines()[-1].split()[-2])
    return simTime


def compare_with_analyt(runInfo, dirDict, tol=1e-4):
    failList = []
    for testStr in sorted(runInfo.keys()):
        newDir = osp.join(dirDict["out"], testStr, "sim_output")
        newDataFile = osp.join(newDir, "output_data.mat")
        config = Config.from_dicts(newDir)

        try:
            newData = sio.loadmat(newDataFile)
        except IOError as exception:
            # If it's an error _other than_ the file not being there
            if exception.errno != errno.ENOENT:
                raise
            print("No simulation data for " + testStr)
            continue
        if "Difn" in testStr:
            t_ref = config["t_ref"]
            L_part = config["psd_len"]["c"][0,0]
            nx_part = config["psd_num"]["c"][0,0]
            t_refPart = L_part**2 / config["c", "D"]
            # Skip first time point: analytical solution fails at t=0.
            t0ind = 2
            r0ind = 1
            tvecA = newData["phi_applied_times"][0][t0ind:] * (t_ref/t_refPart)
            cmat = newData["partTrodecvol0part0_c"]
            cmin, cmax = np.min(cmat), np.max(cmat)
            delC = cmax - cmin
            # Skip center mesh point and first time points:
            # analytical solutions have issues with r=0 and t=0
            cmat = cmat[t0ind:,r0ind:]
            nt = len(tvecA)
            # Skip center mesh point: analytical solution as sin(r)/r
            xvecA = np.linspace(0., 1., nx_part)[1:]
            R, T = np.meshgrid(xvecA, tvecA)
            theta = runInfo[testStr][1](R, T)
            theta = delC*theta + cmin
            for tind in range(nt):
                cvec = cmat[tind,:]
                thetavec = theta[tind,:]
                if np.max(np.abs(thetavec - cvec)) > tol:
                    print(testStr, "Fail from tolerance")
                    failList.append(testStr)
                    break
    return failList


def compare_with_ref(runInfo, dirDict, tol=1e-4):
    timeList_new = []
    timeList_ref = []
    failList = []
    for testStr in sorted(runInfo.keys()):
        testFailed = False
        newDir = osp.join(dirDict["out"], testStr, "sim_output")
        refDir = osp.join(dirDict["refs"], testStr, "sim_output")
        newDataFile = osp.join(newDir, "output_data.mat")
        refDataFile = osp.join(refDir, "output_data.mat")
        timeList_new.append(get_sim_time(newDir))
        timeList_ref.append(get_sim_time(refDir))
        try:
            newData = sio.loadmat(newDataFile)
        except IOError as exception:
            # If it's an error _other than_ the file not being there
            if exception.errno != errno.ENOENT:
                raise
            print("No simulation data for " + testStr)
            continue

        refData = sio.loadmat(refDataFile)
        for varKey in (set(refData.keys()) & set(newData.keys())):
            # TODO -- Consider keeping a list of the variables that fail

            # Ignore certain entries not of numerical output
            if varKey[0:2] == "__":
                continue

            # Compute the difference between the solution and the reference
            try:
                varDataNew = newData[varKey]
                varDataRef = refData[varKey]
                diffMat = np.abs(varDataNew - varDataRef)
            except ValueError:
                print(testStr, "Fail from ValueError")
                testFailed = True
                continue
            except KeyError:
                print(testStr, "Fail from KeyError")
                testFailed = True
                continue

            # #Check absolute and relative error against tol
            if np.max(diffMat) > tol and np.max(diffMat) > tol*np.max(np.abs(varDataRef)):
                print(testStr, "Fail from tolerance")
                print("variable failing:", varKey)
                print("max error:", np.max(diffMat))
                testFailed = True

        if testFailed:
            failList.append(testStr)

    scl = 1.3
    fig, ax = plt.subplots(figsize=(scl*6, scl*4))
    ax.plot(timeList_ref, timeList_new, 'o')
    tmax = max(max(timeList_ref), max(timeList_new))
    ax.plot([0., tmax], [0., tmax], linestyle="--")
    ax.set_xlabel("Reference simulation time [s]")
    ax.set_ylabel("New test simulation time [s]")
    fig.savefig(osp.join(dirDict["plots"], "timeParity.png"),
                bbox_inches="tight")
    plt.close("all")
    return failList


def show_fails(failList):
    for fail in failList:
        print((fail + " differs from the reference outputs!"))


def main(compareDir):
    pflag = True
    dirDict = {}
    # Get the default configs
    dirDict["suite"] = osp.dirname(osp.abspath(__file__))
    dirDict["simOut"] = osp.join(os.getcwd(), "sim_output")
    dirDict["out"] = osp.join(dirDict["suite"], "test_outputs", time.strftime("%Y%m%d_%H%M%S"))
    dirDict["baseConfig"] = osp.join(dirDict["suite"], "baseConfigs")
    dirDict["refs"] = osp.join(dirDict["suite"], "ref_outputs")

    # Dictionary containing info about the tests to run
    # Identifier strings are associated with functions to call and
    # whether to run that particular test.
    n_tests = 19
    runInfo = {'test{:03}'.format(i): getattr(defs, 'test{:03}'.format(i))
               for i in range(1, n_tests+1)}
    runInfoAnalyt = {
        "testAnalytCylDifn": (defs.testAnalytCylDifn, defs.analytCylDifn),
        "testAnalytSphDifn": (defs.testAnalytSphDifn, defs.analytSphDifn),
        }

    if compareDir is None:
        os.makedirs(dirDict["out"])
        dirDict["plots"] = osp.join(dirDict["out"], "plots")
        os.makedirs(dirDict["plots"])
        run_test_sims(runInfo, dirDict, pflag)
        run_test_sims_analyt(runInfoAnalyt, dirDict)
    else:
        dirDict["out"] = compareDir
        dirDict["plots"] = osp.join(dirDict["out"], "plots")
    failList = compare_with_ref(runInfo, dirDict, tol=1e-4)
    failListAnalyt = compare_with_analyt(runInfoAnalyt, dirDict, tol=1e-4)

    # Print a summary of test results
    print("\n")
    print("--------------------------------------------------------------")
    if len(failList) > 0:
        print("Comparison fails:")
        for fail in failList:
            print("\t" + fail + " differs from the reference outputs!")
    else:
        print("All comparison tests passed!")
    if len(failListAnalyt) > 0:
        print("Analytical fails:")
        for fail in failListAnalyt:
            print("\t" + fail + " differs from the reference outputs!")
    else:
        print("All analytical tests passed!")
    print("--------------------------------------------------------------")
