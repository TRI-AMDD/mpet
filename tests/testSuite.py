import os
import os.path as osp
import sys
import errno
import time
import ConfigParser

import numpy as np
import scipy.io as sio
import matplotlib as mpl
mpl.use("TkAgg")
import matplotlib.pyplot as plt
# Plot defaults
axtickfsize = 18
labelfsize = 20
legfsize = labelfsize - 2
txtfsize = labelfsize - 2
lwidth = 3.
markersize = 10
mpl.rcParams['xtick.labelsize'] = axtickfsize
mpl.rcParams['ytick.labelsize'] = axtickfsize
mpl.rcParams['axes.labelsize'] = labelfsize
mpl.rcParams['axes.labelsize'] = labelfsize
mpl.rcParams['font.size'] = txtfsize
mpl.rcParams['legend.fontsize'] = legfsize
mpl.rcParams['lines.linewidth'] = lwidth
mpl.rcParams['lines.markersize'] = markersize
mpl.rcParams['lines.markeredgewidth'] = 0.1
#mpl.rcParams['text.usetex'] = True

mpetdir = osp.join(osp.dirname(osp.abspath(__file__)), "..")
sys.path.append(mpetdir)
import mpetParamsIO as IO
import testDefns as defs
import plot_data as pd

def run_test_sims(runInfo, dirDict, pflag=True):
    for testStr in sorted(runInfo.keys()):
        testDir = osp.join(dirDict["out"], testStr)
        os.makedirs(testDir)
        runInfo[testStr](testDir, dirDict, pflag)
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
        dD_s, ndD_s = IO.readDicts(osp.join(newDir, "input_dict_system"))
        tmp = IO.readDicts(osp.join(newDir, "input_dict_c"))
        dD_e = {}
        ndD_e = {}
        dD_e["c"], ndD_e["c"] = tmp
        try:
            newData = sio.loadmat(newDataFile)
        except IOError as exception:
            # If it's an error _other than_ the file not being there
            if exception.errno != errno.ENOENT:
                raise
            print "No simulation data for " + testStr
            continue
        if "Difn" in testStr:
            t_ref = dD_s["t_ref"]
            L_part = dD_s["psd_len"]["c"][0, 0]
            nx_part = ndD_s["psd_num"]["c"][0, 0]
            t_refPart = L_part**2 / dD_e["c"]["D"]
            # Skip first time point: analytical solution fails at t=0.
            t0ind = 2
            r0ind = 1
            tvecA = newData["phi_applied_times"][0][t0ind:] * (t_ref/t_refPart)
            cmat = newData["partTrodecvol0part0_c"]
            cmin, cmax = np.min(cmat), np.max(cmat)
            delC = cmax - cmin
            # Skip center mesh point and first time points:
            # analytical solutions have issues with r=0 and t=0
            cmat = cmat[t0ind:, r0ind:]
            nt = len(tvecA)
            # Skip center mesh point: analytical solution as sin(r)/r
            xvecA = np.linspace(0., 1., nx_part)[1:]
            R, T = np.meshgrid(xvecA, tvecA)
            theta = runInfo[testStr][1](R, T)
            theta = delC*theta + cmin
            for tind in range(nt):
                cvec = cmat[tind, :]
                thetavec = theta[tind, :]
                if np.max(np.abs(thetavec - cvec)) > tol:
                    print testStr, "Fail from tolerance"
                    failList.append(testStr)
                    break
    return failList

def compare_with_ref(runInfo, dirDict, tol=1e-4):
    timeList_new = []
    timeList_ref = []
    failList = []
    for testStr in sorted(runInfo.keys()):
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
            print "No simulation data for " + testStr
            continue
        refData = sio.loadmat(refDataFile)
        for varKey in refData.keys():
            # If this test has already failed
            # TODO -- Consider keeping a list of the variables that fail
            if testStr in failList:
                break
            # Ignore certain entries not of numerical output
            if varKey[0:2] == "__":
                continue
            try:
                varDataNew = newData[varKey]
                varDataRef = refData[varKey]
                diffMat = np.abs(varDataNew - varDataRef)
            except ValueError:
                print testStr, "Fail from ValueError"
                failList.append(testStr)
                continue
            except KeyError:
                print testStr, "Fail from KeyError"
                failList.append(testStr)
                continue
            # TODO -- What is the right way to compare here?
            absTol = tol*(np.max(varDataNew) - np.min(varDataNew))
            if np.max(diffMat) > absTol:
                print testStr, "Fail from tolerance"
                print "variable failing:", varKey
                print "max error:", np.max(diffMat)
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
        print (fail + " differs from the reference outputs!")

def main(compareDir):
    pflag = True
    dirDict = {}
    # Get the default configs
    dirDict["suite"] = osp.dirname(osp.abspath(__file__))
    dirDict["simOut"] = osp.join(os.getcwd(), "sim_output")
    dirDict["out"] = osp.join(dirDict["suite"],
                      time.strftime("%Y%m%d_%H%M%S", time.localtime()))
    dirDict["baseConfig"] = osp.join(dirDict["suite"], "baseConfigs")
    dirDict["refs"] = osp.join(dirDict["suite"], "ref_outputs")

    # Dictionary containing info about the tests to run
    # Identifier strings are associated with functions to call and
    # whether to run that particular test.
    runInfo = {
            "test001" : defs.test001,
            "test002" : defs.test002,
            "test003" : defs.test003,
            "test004" : defs.test004,
            "test005" : defs.test005,
            "test006" : defs.test006,
            "test007" : defs.test007,
            "test008" : defs.test008,
            "test009" : defs.test009,
            "test010" : defs.test010,
            "test011" : defs.test011,
            "test012" : defs.test012,
            "test013" : defs.test013,
            "test014" : defs.test014,
            "test015" : defs.test015,
            "test016" : defs.test016,
            "test017" : defs.test017,
            "test018" : defs.test018,
            }
    runInfoAnalyt = {
            "testAnalytCylDifn" : (defs.testAnalytCylDifn, defs.analytCylDifn),
            "testAnalytSphDifn" : (defs.testAnalytSphDifn, defs.analytSphDifn),
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
    failList = compare_with_ref(runInfo, dirDict, tol=1e-3)
    failListAnalyt = compare_with_analyt(runInfoAnalyt, dirDict, tol=1e-4)

    if len(failList) > 0:
        print "Comparison fails:"
        show_fails(failList)
    else:
        print "All comparison tests passed!"
    if len(failListAnalyt) > 0:
        print "Analytical fails:"
        show_fails(failListAnalyt)
    else:
        print "All analytical tests passed!"

if __name__ == "__main__":
    if len(sys.argv) > 1:
        compareDir = osp.join(os.getcwd(), sys.argv[1])
    else:
        compareDir = None
    main(compareDir)
