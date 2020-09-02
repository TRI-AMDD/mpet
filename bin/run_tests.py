import os.path as osp
from os import getcwd, makedirs,rmdir
from sys import argv
from tests.test_suite import run_test_sims, run_test_sims_analyt
import tests.test_defs as defs
from shutil import rmtree

def run(dirName):
    pflag = False
    dirDict = {}
    # Get the default configs
    dirDict["suite"] = osp.join(osp.dirname(osp.abspath(__file__)),"tests")
    dirDict["simOut"] = osp.join(getcwd(), "sim_output")
    dirDict["out"] = osp.join(dirDict["suite"], "test_outputs", dirName)
    dirDict["baseConfig"] = osp.join(dirDict["suite"], "baseConfigs")

    # Dictionary containing info about the tests to run
    # Identifier strings are associated with functions to call and
    # whether to run that particular test.
    n_tests = defs.getNumberOfTests()
    runInfo = {'test{:03}'.format(i): getattr(defs, 'test{:03}'.format(i))
               for i in range(1, n_tests+1)}
    runInfoAnalyt = {
        "testAnalytCylDifn": (defs.testAnalytCylDifn, defs.analytCylDifn),
        "testAnalytSphDifn": (defs.testAnalytSphDifn, defs.analytSphDifn),
        }

    if osp.exists(dirDict["out"]):
      rmtree(dirDict["out"])
    makedirs(dirDict["out"])
    dirDict["plots"] = osp.join(dirDict["out"], "plots")
    makedirs(dirDict["plots"])
    run_test_sims(runInfo, dirDict, pflag)
    run_test_sims_analyt(runInfoAnalyt, dirDict)

if __name__ == "__main__":
  run(argv[1])
