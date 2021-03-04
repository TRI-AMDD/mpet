import os.path as osp
import pytest


def get_sim_time(simDir):
    with open(osp.join(simDir, "run_info.txt")) as fi:
        simTime = float(fi.readlines()[-1].split()[-2])
    return simTime


def test_compare_timings(Dirs, tol):
    refDir, testDir = Dirs
    newDir = osp.join(testDir, "sim_output")
    refDir = osp.join(refDir, "sim_output")
    time_new = get_sim_time(newDir)
    time_ref = get_sim_time(refDir)
    assert time_new < time_ref * pytest.slowdown_tolerance, "too high slowdown with reference"
