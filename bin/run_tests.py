import os.path as osp
import time
from os import getcwd, makedirs,rmdir
from sys import argv
from tests.test_suite import run_test_sims, run_test_sims_analyt
import tests.test_defs as defs
from shutil import rmtree
import argparse

def run(test_outputs, testDir, tests=None):
    pflag = False
    dirDict = {}
    # Get the default configs
    dirDict["suite"] = testDir
    dirDict["simOut"] = osp.join(getcwd(), "sim_output")
    dirDict["out"] = test_outputs
    dirDict["baseConfig"] = osp.join(dirDict["suite"], "baseConfigs")

    # Dictionary containing info about the tests to run
    # Identifier strings are associated with functions to call and
    # whether to run that particular test.
    n_tests = defs.getNumberOfTests()
    if tests:
      runInfo = { test : getattr(defs, test) for test in tests }
    else:
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
    if not tests:
      run_test_sims_analyt(runInfoAnalyt, dirDict)

def main():
    parser = argparse.ArgumentParser(description='Run test suite.')
    parser.add_argument('--output_dir', metavar='o', type=str,
        default=osp.join(osp.dirname(osp.abspath(__file__)),"../tests",
          "test_outputs", time.strftime("%Y%m%d_%H%M%S") ),
                      help='where to put the simulations, default is tests/tests_outputs/date')
    parser.add_argument('--test_dir', metavar='t', type=str, 
        default=osp.join(osp.dirname(osp.abspath(__file__)),"../tests"),
                      help='where are the tests located?')
    parser.add_argument('tests', nargs='*', default=[], help='which tests do I run?')
    args = parser.parse_args()

    run(args.output_dir,args.test_dir, args.tests)

    return(args)

if __name__ == "__main__":
    main()
