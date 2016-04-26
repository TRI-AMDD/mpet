import os
import os.path as osp
import sys
import errno
import time
import ConfigParser

import numpy as np
import scipy.io as sio

mpetdir = osp.join(os.environ["HOME"], "docs", "bazantgroup", "mpet")
sys.path.append(mpetdir)
import mpetParamsIO as IO
import testDefns as defs

def run_test_sims(runInfo, dirDict):
    for testStr in sorted(runInfo.keys()):
        testDir = osp.join(dirDict["out"], testStr)
        os.makedirs(testDir)
        runInfo[testStr](testDir, dirDict)
    # Remove the history directory that mpet creates.
    try:
        os.rmdir(osp.join(dirDict["suite"], "history"))
    except OSError as exception:
        if exception.errno != errno.ENOENT:
            raise
    return

def compare_with_ref(runInfo, dirDict, tol=1e-4):
    failList = []
    for testStr in sorted(runInfo.keys()):
        newDataFile = osp.join(
            dirDict["out"], testStr, "sim_output", "output_data.mat")
        refDataFile = osp.join(
            dirDict["refs"], testStr, "sim_output", "output_data.mat")
        try:
            newData = sio.loadmat(newDataFile)
        except IOError as exception:
            # If it's an error _other than_ the file not being there
            if exception.errno != errno.ENOENT:
                raise
            print "No simulation data for " + testStr
            continue
        refData = sio.loadmat(refDataFile)
        for varKey in newData.keys():
            # If this test has already failed
            # TODO -- Consider keeping a list of the variables that fail
            if testStr in failList:
                continue
            # Ignore certain entries not of numerical output
            if varKey[0:2] == "__":
                continue
            varDataNew = newData[varKey]
            varDataRef = refData[varKey]
            try:
                diffMat = np.abs(varDataNew - varDataRef)
            except ValueError:
                failList.append(testStr)
                continue
            # TODO -- What is the right way to compare here?
            absTol = tol*(np.max(varDataNew) - np.min(varDataNew))
            if np.max(diffMat) > absTol:
                failList.append(testStr)
    return failList

def show_fails(failList):
    for fail in failList:
        print (fail + " differs from the reference outputs!")
    return

def main(compareDir):
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
        run_test_sims(runInfo, dirDict)
    else:
        dirDict["out"] = compareDir
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
