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

mpetdir = osp.join(os.environ["HOME"], "docs", "bazantgroup", "mpet")
sys.path.append(mpetdir)
import mpetParamsIO as IO
import testDefns as defs

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
    return

def get_sim_time(simDir):
    with open(osp.join(simDir, "run_info.txt")) as fi:
        simTime = float(fi.readlines()[-1].split()[-2])
    return simTime

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
            varDataNew = newData[varKey]
            varDataRef = refData[varKey]
            try:
                diffMat = np.abs(varDataNew - varDataRef)
            except ValueError:
                print testStr, "Fail from ValueError"
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
    return

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

    if compareDir is None:
        os.makedirs(dirDict["out"])
        dirDict["plots"] = osp.join(dirDict["out"], "plots")
        os.makedirs(dirDict["plots"])
        run_test_sims(runInfo, dirDict, pflag)
    else:
        dirDict["out"] = compareDir
        dirDict["plots"] = osp.join(dirDict["out"], "plots")
    failList = compare_with_ref(runInfo, dirDict, tol=1e-3)

    if len(failList) > 0:
        show_fails(failList)
    else:
        print "All tests passed!"

    return

if __name__ == "__main__":
    if len(sys.argv) > 1:
        compareDir = osp.join(os.getcwd(), sys.argv[1])
    else:
        compareDir = None
    main(compareDir)
