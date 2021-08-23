import os.path as osp
import shutil

import configparser
import numpy as np
import scipy.special as spcl

import tests.compare_plots as cmpr
import mpet.main as main
import re


def corePlots(testDir, dirDict):
    cmpr.vt(testDir, dirDict)
    cmpr.curr(testDir, dirDict)


def elytePlots(testDir, dirDict):
    cmpr.elytecons(testDir, dirDict)
    cmpr.elytecf(testDir, dirDict)
    cmpr.elytepf(testDir, dirDict)
    cmpr.elyteif(testDir, dirDict)
    cmpr.elytedivif(testDir, dirDict)


def electrodePlots(testDir, dirDict, trode):
    cmpr.soc(testDir, dirDict, "c")
    cmpr.cbarLine(testDir, dirDict, "c")


def getNumberOfTests():
    n = 0
    for key in globals().keys():
        if re.match(r'test[0-9][0-9][0-9]$', key):
            n = n + 1
    return n


def get_config(inFile):
    P = configparser.RawConfigParser()
    P.optionxform = str
    P.read(inFile)
    return P


def write_config_file(P, filename="input_params.cfg"):
    with open(filename, "w") as fo:
        P.write(fo)


def testAnalytSphDifn(testDir, dirDict):
    """ Analytical test for diffusion in a sphere """
    shutil.copy(osp.join(dirDict["baseConfig"], "params_system.cfg"), testDir)
    shutil.copy(osp.join(dirDict["baseConfig"], "params_c.cfg"), testDir)
    psys = osp.join(testDir, "params_system.cfg")
    ptrode = osp.join(testDir, "params_c.cfg")

    P_s = get_config(psys)
    P_s.set("Sim Params", "profileType", "CV")
    P_s.set("Sim Params", "Vset", "0.10")
    P_s.set("Sim Params", "tend", "1e+6")
    P_s.set("Sim Params", "tsteps", "1800")
    P_s.set("Sim Params", "relTol", "1e-7")
    P_s.set("Sim Params", "tramp", "1e-6")
    P_s.set("Particles", "mean_c", "100e-9")
    P_s.set("Particles", "mean_c", "100e-9")
    P_s.set("Particles", "cs0_c", "0.5")
    write_config_file(P_s, psys)
    P = get_config(ptrode)
    P.set("Particles", "type", "diffn")
    P.set("Particles", "discretization", "3e-10")
    P.set("Particles", "shape", "sphere")
    P.set("Material", "muRfunc", "testIS_ss")
    P.set("Material", "D", "1e-20")
    P.set("Reactions", "rxnType", "BV_raw")
    P.set("Reactions", "k0", "1e+1")
    write_config_file(P, ptrode)
    main.main(psys, keepArchive=False)
    shutil.move(dirDict["simOut"], testDir)


def testAnalytCylDifn(testDir, dirDict):
    """ Analytical test for diffusion in a cylinder """
    shutil.copy(osp.join(dirDict["baseConfig"], "params_system.cfg"), testDir)
    shutil.copy(osp.join(dirDict["baseConfig"], "params_c.cfg"), testDir)
    psys = osp.join(testDir, "params_system.cfg")
    ptrode = osp.join(testDir, "params_c.cfg")
    P_s = get_config(psys)
    P_s.set("Sim Params", "profileType", "CV")
    P_s.set("Sim Params", "Vset", "0.10")
    P_s.set("Sim Params", "tend", "1e+6")
    P_s.set("Sim Params", "tsteps", "1800")
    P_s.set("Sim Params", "relTol", "1e-7")
    P_s.set("Sim Params", "tramp", "1e-6")
    P_s.set("Particles", "mean_c", "100e-9")
    P_s.set("Particles", "mean_c", "100e-9")
    P_s.set("Particles", "cs0_c", "0.5")
    write_config_file(P_s, psys)
    P = get_config(ptrode)
    P.set("Particles", "type", "diffn")
    P.set("Particles", "discretization", "3e-10")
    P.set("Particles", "shape", "cylinder")
    P.set("Material", "muRfunc", "testIS_ss")
    P.set("Material", "D", "1e-20")
    P.set("Reactions", "rxnType", "BV_raw")
    P.set("Reactions", "k0", "1e+1")
    write_config_file(P, ptrode)
    main.main(psys, keepArchive=False)
    shutil.move(dirDict["simOut"], testDir)


def analytSphDifn(R, T):
    """ Analytical solution for diffusion in a sphere """
    n = 30
    lmbdavec = np.pi * np.arange(1, n + 1)
    theta = 0 * R
    for i, lmbda in enumerate(lmbdavec):
        theta += ((2. / lmbda ** 2) * np.sin(lmbda * R) / R
                  * (np.sin(lmbda) - lmbda * np.cos(lmbda))
                  * np.exp(-lmbda ** 2 * T))
    return theta


def analytCylDifn(R, T):
    """ Analytical solution for diffusion in a cylinder """
    n = 30
    lmbdavec = spcl.jn_zeros(0, n)
    theta = 0 * R
    for i, lmbda in enumerate(lmbdavec):
        theta += ((2. / lmbda) * spcl.j0(lmbda * R) / spcl.j1(lmbda)
                  * np.exp(-lmbda ** 2 * T))
    return theta
