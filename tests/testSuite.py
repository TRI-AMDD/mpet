import os
import os.path as osp
import sys
import errno
import shutil
import time
import ConfigParser

import numpy as np
import scipy.io as sio

mpetdir = osp.join(os.environ["HOME"], "docs", "bazantgroup", "mpet")
sys.path.append(mpetdir)
import mpet
import mpetParamsIO

def test001(IO, testDir, baseConfigDir, simOutDir, run=True):
    """ LFP ACR C3 """
    shutil.copy(osp.join(baseConfigDir, "params_system.cfg"), testDir)
    shutil.copy(osp.join(baseConfigDir, "params_c.cfg"), testDir)
    psys = osp.join(testDir, "params_system.cfg")
    ptrode = osp.join(testDir, "params_c.cfg")
    P_s = mpetParamsIO.getConfig(psys)
    P_s.set("Sim Params", "Crate", "1e-2")
    IO.writeConfigFile(P_s, psys)
    P = mpetParamsIO.getConfig(ptrode)
    P.set("Particles", "type", "ACR")
    P.set("Material", "muRfunc", "LiFePO4")
    IO.writeConfigFile(P, ptrode)
    if run:
        mpet.main(psys, keepArchive=False)
        shutil.move(simOutDir, osp.join(testDir))
    return

def test002(IO, testDir, baseConfigDir, simOutDir, run=True):
    """ LFP CHR cylinder """
    shutil.copy(osp.join(baseConfigDir, "params_system.cfg"), testDir)
    shutil.copy(osp.join(baseConfigDir, "params_c.cfg"), testDir)
    psys = osp.join(testDir, "params_system.cfg")
    ptrode = osp.join(testDir, "params_c.cfg")
    P_s = mpetParamsIO.getConfig(psys)
    P_s.set("Sim Params", "Crate", "1e-2")
    IO.writeConfigFile(P_s, psys)
    P = mpetParamsIO.getConfig(ptrode)
    P.set("Particles", "type", "CHR")
    P.set("Particles", "shape", "cylinder")
    P.set("Material", "muRfunc", "LiFePO4")
    P.set("Material", "dgammadc", "5e-30")
    IO.writeConfigFile(P, ptrode)
    if run:
        mpet.main(psys, keepArchive=False)
        shutil.move(simOutDir, osp.join(testDir))
    return

def test003(IO, testDir, baseConfigDir, simOutDir, run=True):
    """ LFP CHR sphere """
    shutil.copy(osp.join(baseConfigDir, "params_system.cfg"), testDir)
    shutil.copy(osp.join(baseConfigDir, "params_c.cfg"), testDir)
    psys = osp.join(testDir, "params_system.cfg")
    ptrode = osp.join(testDir, "params_c.cfg")
    P_s = mpetParamsIO.getConfig(psys)
    P_s.set("Sim Params", "Crate", "1e-2")
    IO.writeConfigFile(P_s, psys)
    P = mpetParamsIO.getConfig(ptrode)
    P.set("Particles", "type", "CHR")
    P.set("Particles", "shape", "sphere")
    P.set("Material", "muRfunc", "LiFePO4")
    P.set("Material", "dgammadc", "-2e-30")
    IO.writeConfigFile(P, ptrode)
    if run:
        mpet.main(psys, keepArchive=False)
        shutil.move(simOutDir, osp.join(testDir))
    return

def test004(IO, testDir, baseConfigDir, simOutDir, run=True):
    """ LFP CHR sphere with noise  """
    shutil.copy(osp.join(baseConfigDir, "params_system.cfg"), testDir)
    shutil.copy(osp.join(baseConfigDir, "params_c.cfg"), testDir)
    psys = osp.join(testDir, "params_system.cfg")
    ptrode = osp.join(testDir, "params_c.cfg")
    P_s = mpetParamsIO.getConfig(psys)
    P_s.set("Sim Params", "Crate", "1e-2")
    IO.writeConfigFile(P_s, psys)
    P = mpetParamsIO.getConfig(ptrode)
    P.set("Particles", "type", "CHR")
    P.set("Particles", "shape", "sphere")
    P.set("Material", "muRfunc", "LiFePO4")
    P.set("Material", "noise", "true")
    P.set("Material", "dgammadc", "-2e-30")
    IO.writeConfigFile(P, ptrode)
    if run:
        mpet.main(psys, keepArchive=False)
        shutil.move(simOutDir, osp.join(testDir))
    return

def test005(IO, testDir, baseConfigDir, simOutDir, run=True):
    """ LFP homog """
    shutil.copy(osp.join(baseConfigDir, "params_system.cfg"), testDir)
    shutil.copy(osp.join(baseConfigDir, "params_c.cfg"), testDir)
    psys = osp.join(testDir, "params_system.cfg")
    ptrode = osp.join(testDir, "params_c.cfg")
    P_s = mpetParamsIO.getConfig(psys)
    P_s.set("Sim Params", "Crate", "1e-2")
    P_s.set("Sim Params", "Nvol_c", "2")
    P_s.set("Sim Params", "Nvol_s", "2")
    IO.writeConfigFile(P_s, psys)
    P = mpetParamsIO.getConfig(ptrode)
    P.set("Particles", "type", "homog")
    P.set("Material", "muRfunc", "LiFePO4")
    IO.writeConfigFile(P, ptrode)
    if run:
        mpet.main(psys, keepArchive=False)
        shutil.move(simOutDir, osp.join(testDir))
    return

def test006(IO, testDir, baseConfigDir, simOutDir, run=True):
    """ LFP homog with logPad, Vmin """
    shutil.copy(osp.join(baseConfigDir, "params_system.cfg"), testDir)
    shutil.copy(osp.join(baseConfigDir, "params_c.cfg"), testDir)
    psys = osp.join(testDir, "params_system.cfg")
    ptrode = osp.join(testDir, "params_c.cfg")
    P_s = mpetParamsIO.getConfig(psys)
    P_s.set("Sim Params", "Crate", "1e-2")
    P_s.set("Sim Params", "Vmin", "2.8e-0")
    P_s.set("Particles", "cs0_c", "0.5")
    IO.writeConfigFile(P_s, psys)
    P = mpetParamsIO.getConfig(ptrode)
    P.set("Particles", "type", "homog")
    P.set("Material", "muRfunc", "LiFePO4")
    IO.writeConfigFile(P, ptrode)
    if run:
        mpet.main(psys, keepArchive=False)
        shutil.move(simOutDir, osp.join(testDir))
    return

def test007(IO, testDir, baseConfigDir, simOutDir, run=True):
    """ LFP homog with logPad, Vmax """
    shutil.copy(osp.join(baseConfigDir, "params_system.cfg"), testDir)
    shutil.copy(osp.join(baseConfigDir, "params_c.cfg"), testDir)
    psys = osp.join(testDir, "params_system.cfg")
    ptrode = osp.join(testDir, "params_c.cfg")
    P_s = mpetParamsIO.getConfig(psys)
    P_s.set("Sim Params", "Crate", "-1e-2")
    P_s.set("Sim Params", "Vmax", "4.0e-0")
    P_s.set("Particles", "cs0_c", "0.5")
    IO.writeConfigFile(P_s, psys)
    P = mpetParamsIO.getConfig(ptrode)
    P.set("Particles", "type", "homog")
    P.set("Material", "muRfunc", "LiFePO4")
    IO.writeConfigFile(P, ptrode)
    if run:
        mpet.main(psys, keepArchive=False)
        shutil.move(simOutDir, osp.join(testDir))
    return

def test008(IO, testDir, baseConfigDir, simOutDir, run=True):
    """ LFP homog_sdn """
    shutil.copy(osp.join(baseConfigDir, "params_system.cfg"), testDir)
    shutil.copy(osp.join(baseConfigDir, "params_c.cfg"), testDir)
    psys = osp.join(testDir, "params_system.cfg")
    ptrode = osp.join(testDir, "params_c.cfg")
    P_s = mpetParamsIO.getConfig(psys)
    P_s.set("Sim Params", "Crate", "1e-2")
    P_s.set("Sim Params", "Nvol_c", "3")
    P_s.set("Sim Params", "Npart_c", "3")
    P_s.set("Particles", "stddev_c", "25e-9")
    P_s.set("Sim Params", "capFrac", "0.6")
    IO.writeConfigFile(P_s, psys)
    P = mpetParamsIO.getConfig(ptrode)
    P.set("Particles", "type", "homog_sdn")
    P.set("Material", "muRfunc", "LiFePO4")
    IO.writeConfigFile(P, ptrode)
    if run:
        mpet.main(psys, keepArchive=False)
        shutil.move(simOutDir, osp.join(testDir))
    return

def test009(IO, testDir, baseConfigDir, simOutDir, run=True):
    """ Graphite-2param homog """
    shutil.copy(osp.join(baseConfigDir, "params_system.cfg"), testDir)
    shutil.copy(osp.join(baseConfigDir, "params_c.cfg"), testDir)
    psys = osp.join(testDir, "params_system.cfg")
    ptrode = osp.join(testDir, "params_c.cfg")
    P_s = mpetParamsIO.getConfig(psys)
    P_s.set("Sim Params", "Crate", "1e-2")
    IO.writeConfigFile(P_s, psys)
    P = mpetParamsIO.getConfig(ptrode)
    P.set("Particles", "type", "homog2")
    P.set("Material", "muRfunc", "LiC6")
    IO.writeConfigFile(P, ptrode)
    if run:
        mpet.main(psys, keepArchive=False)
        shutil.move(simOutDir, osp.join(testDir))
    return

def test010(IO, testDir, baseConfigDir, simOutDir, run=True):
    """ Graphite-2param CHR cylinder """
    shutil.copy(osp.join(baseConfigDir, "params_system.cfg"), testDir)
    shutil.copy(osp.join(baseConfigDir, "params_c.cfg"), testDir)
    psys = osp.join(testDir, "params_system.cfg")
    ptrode = osp.join(testDir, "params_c.cfg")
    P_s = mpetParamsIO.getConfig(psys)
    P_s.set("Sim Params", "Crate", "1e-2")
    IO.writeConfigFile(P_s, psys)
    P = mpetParamsIO.getConfig(ptrode)
    P.set("Particles", "type", "CHR2")
    P.set("Particles", "shape", "cylinder")
    P.set("Material", "muRfunc", "LiC6")
    IO.writeConfigFile(P, ptrode)
    if run:
        mpet.main(psys, keepArchive=False)
        shutil.move(simOutDir, osp.join(testDir))
    return

def test011(IO, testDir, baseConfigDir, simOutDir, run=True):
    """ Graphite-2param CHR sphere """
    shutil.copy(osp.join(baseConfigDir, "params_system.cfg"), testDir)
    shutil.copy(osp.join(baseConfigDir, "params_c.cfg"), testDir)
    psys = osp.join(testDir, "params_system.cfg")
    ptrode = osp.join(testDir, "params_c.cfg")
    P_s = mpetParamsIO.getConfig(psys)
    P_s.set("Sim Params", "Crate", "1e-2")
    P_s.set("Sim Params", "relTol", "1e-7")
    P_s.set("Sim Params", "absTol", "1e-7")
    IO.writeConfigFile(P_s, psys)
    P = mpetParamsIO.getConfig(ptrode)
    P.set("Particles", "type", "CHR2")
    P.set("Particles", "shape", "sphere")
    P.set("Material", "muRfunc", "LiC6")
    IO.writeConfigFile(P, ptrode)
    if run:
        mpet.main(psys, keepArchive=False)
        shutil.move(simOutDir, osp.join(testDir))
    return

def test012(IO, testDir, baseConfigDir, simOutDir, run=True):
    """ Solid solution, diffn sphere, homog, LiC6_coke_ss2, LiMn2O4_ss2
    BV_mod01, BV_mod02
    cathode + anode
    """
    shutil.copy(osp.join(baseConfigDir, "params_system.cfg"), testDir)
    shutil.copy(osp.join(baseConfigDir, "params_c.cfg"), testDir)
    shutil.copy(osp.join(baseConfigDir, "params_a.cfg"), testDir)
    psys = osp.join(testDir, "params_system.cfg")
    ptrodec = osp.join(testDir, "params_c.cfg")
    ptrodea = osp.join(testDir, "params_a.cfg")
    P_s = mpetParamsIO.getConfig(psys)
    P_s.set("Sim Params", "Crate", "1e-2")
    P_s.set("Sim Params", "capFrac", "0.67")
    P_s.set("Sim Params", "Nvol_c", "2")
    P_s.set("Sim Params", "Nvol_a", "2")
    P_s.set("Sim Params", "Vmin", "2e0")
    P_s.set("Particles", "cs0_c", "0.2")
    P_s.set("Particles", "cs0_a", "0.495")
    P_s.set("Electrolyte", "elyteModelType", "SM")
    IO.writeConfigFile(P_s, psys)
    P = mpetParamsIO.getConfig(ptrodec)
    P.set("Particles", "type", "homog")
    P.set("Particles", "shape", "cylinder")
    P.set("Material", "muRfunc", "LiMn2O4_ss2")
    P.set("Reactions", "rxnType", "BV_mod01")
    IO.writeConfigFile(P, ptrodec)
    P = mpetParamsIO.getConfig(ptrodea)
    P.set("Particles", "shape", "sphere")
    P.set("Material", "muRfunc", "LiC6_coke_ss2")
    P.set("Reactions", "rxnType", "BV_mod02")
    IO.writeConfigFile(P, ptrodea)
    if run:
        mpet.main(psys, keepArchive=False)
        shutil.move(simOutDir, osp.join(testDir))
    return

def test013(IO, testDir, baseConfigDir, simOutDir, run=True):
    """ Solid solution, diffn cylinder, homog, testIS_ss, LiMn2O4_ss2
    Marcus, BV_raw
    cathode + separator + anode
    """
    shutil.copy(osp.join(baseConfigDir, "params_system.cfg"), testDir)
    shutil.copy(osp.join(baseConfigDir, "params_c.cfg"), testDir)
    shutil.copy(osp.join(baseConfigDir, "params_a.cfg"), testDir)
    psys = osp.join(testDir, "params_system.cfg")
    ptrodec = osp.join(testDir, "params_c.cfg")
    ptrodea = osp.join(testDir, "params_a.cfg")
    P_s = mpetParamsIO.getConfig(psys)
    P_s.set("Sim Params", "profileType", "CV")
    P_s.set("Sim Params", "Vset", "3.8")
    P_s.set("Sim Params", "Nvol_c", "2")
    P_s.set("Sim Params", "Nvol_s", "2")
    P_s.set("Sim Params", "Nvol_a", "2")
    P_s.set("Particles", "cs0_c", "0.2")
    P_s.set("Particles", "cs0_a", "0.95")
    P_s.set("Electrolyte", "elyteModelType", "dilute")
    IO.writeConfigFile(P_s, psys)
    P = mpetParamsIO.getConfig(ptrodec)
    P.set("Particles", "type", "homog")
    P.set("Particles", "shape", "sphere")
    P.set("Material", "muRfunc", "LiMn2O4_ss2")
    P.set("Reactions", "rxnType", "Marcus")
    IO.writeConfigFile(P, ptrodec)
    P = mpetParamsIO.getConfig(ptrodea)
    P.set("Particles", "shape", "cylinder")
    P.set("Material", "muRfunc", "testIS_ss")
    P.set("Reactions", "rxnType", "BV_raw")
    IO.writeConfigFile(P, ptrodea)
    if run:
        mpet.main(psys, keepArchive=False)
        shutil.move(simOutDir, osp.join(testDir))
    return

def test014(IO, testDir, baseConfigDir, simOutDir, run=True):
    """ LFP homog with CCsegments, MHC, Rser """
    shutil.copy(osp.join(baseConfigDir, "params_system.cfg"), testDir)
    shutil.copy(osp.join(baseConfigDir, "params_c.cfg"), testDir)
    psys = osp.join(testDir, "params_system.cfg")
    ptrode = osp.join(testDir, "params_c.cfg")
    P_s = mpetParamsIO.getConfig(psys)
    P_s.set("Sim Params", "profileType", "CCsegments")
    P_s.set("Sim Params", "segments",
            "[(1., 25), (-2., 10), (0., 30)]")
    P_s.set("Sim Params", "Rser", "1e-3")
    IO.writeConfigFile(P_s, psys)
    P = mpetParamsIO.getConfig(ptrode)
    P.set("Particles", "type", "homog")
    P.set("Material", "muRfunc", "LiFePO4")
    P.set("Reactions", "rxnType", "MHC")
    IO.writeConfigFile(P, ptrode)
    if run:
        mpet.main(psys, keepArchive=False)
        shutil.move(simOutDir, osp.join(testDir))
    return

def test015(IO, testDir, baseConfigDir, simOutDir, run=True):
    """ testRS homog with CVsegments, bulkCond, partCond """
    shutil.copy(osp.join(baseConfigDir, "params_system.cfg"), testDir)
    shutil.copy(osp.join(baseConfigDir, "params_c.cfg"), testDir)
    psys = osp.join(testDir, "params_system.cfg")
    ptrode = osp.join(testDir, "params_c.cfg")
    P_s = mpetParamsIO.getConfig(psys)
    P_s.set("Sim Params", "profileType", "CVsegments")
    P_s.set("Sim Params", "Npart_c", "5")
    P_s.set("Sim Params", "Nvol_c", "5")
    P_s.set("Sim Params", "segments",
            "[(-0.3, 25), (0., 10), (0.3, 30)]")
    P_s.set("Conductivity", "simBulkCond_c", "true")
    P_s.set("Conductivity", "mcond_c", "5e-2")
    P_s.set("Conductivity", "simPartCond_c", "true")
    P_s.set("Conductivity", "G_mean_c", "3e-14")
    IO.writeConfigFile(P_s, psys)
    P = mpetParamsIO.getConfig(ptrode)
    P.set("Particles", "type", "homog")
    P.set("Material", "muRfunc", "testRS")
    P.set("Reactions", "rxnType", "BV_raw")
    IO.writeConfigFile(P, ptrode)
    if run:
        mpet.main(psys, keepArchive=False)
        shutil.move(simOutDir, osp.join(testDir))
    return

def test016(IO, testDir, baseConfigDir, simOutDir, run=True):
    """ test CC continuation """
    shutil.copy(osp.join(baseConfigDir, "params_system.cfg"), testDir)
    shutil.copy(osp.join(baseConfigDir, "params_c.cfg"), testDir)
    psys = osp.join(testDir, "params_system.cfg")
    ptrode = osp.join(testDir, "params_c.cfg")
    P_s = mpetParamsIO.getConfig(psys)
    P_s.set("Sim Params", "Crate", "1e-2")
    P_s.set("Sim Params", "Nvol_c", "3")
    P_s.set("Sim Params", "Npart_c", "3")
    P_s.set("Sim Params", "capFrac", "0.34")
    test008dir = str(osp.join(testDir, "..", "test008", "sim_output"))
    P_s.set("Sim Params", "prevDir", test008dir)
    IO.writeConfigFile(P_s, psys)
    P = mpetParamsIO.getConfig(ptrode)
    P.set("Particles", "type", "homog_sdn")
    P.set("Material", "muRfunc", "LiFePO4")
    IO.writeConfigFile(P, ptrode)
    if run:
        mpet.main(psys, keepArchive=False)
        shutil.move(simOutDir, osp.join(testDir))
    return

def test017(IO, testDir, baseConfigDir, simOutDir, run=True):
    """ test CV continuation """
    shutil.copy(osp.join(baseConfigDir, "params_system.cfg"), testDir)
    shutil.copy(osp.join(baseConfigDir, "params_c.cfg"), testDir)
    psys = osp.join(testDir, "params_system.cfg")
    ptrode = osp.join(testDir, "params_c.cfg")
    P_s = mpetParamsIO.getConfig(psys)
    P_s.set("Sim Params", "profileType", "CV")
    P_s.set("Sim Params", "Vset", "3.45")
    P_s.set("Sim Params", "tend", "3e3")
    P_s.set("Sim Params", "Nvol_c", "3")
    P_s.set("Sim Params", "Npart_c", "3")
    test008dir = str(osp.join(testDir, "..", "test008", "sim_output"))
    P_s.set("Sim Params", "prevDir", test008dir)
    IO.writeConfigFile(P_s, psys)
    P = mpetParamsIO.getConfig(ptrode)
    P.set("Particles", "type", "homog_sdn")
    P.set("Material", "muRfunc", "LiFePO4")
    IO.writeConfigFile(P, ptrode)
    if run:
        mpet.main(psys, keepArchive=False)
        shutil.move(simOutDir, osp.join(testDir))
    return

def run_test_sims(IO, runInfo, dirDict):
    for testStr in sorted(runInfo.keys()):
        testDir = osp.join(dirDict["out"], testStr)
        os.makedirs(testDir)
        runInfo[testStr](
            IO, testDir, dirDict["baseConfig"], dirDict["simOut"],
            run=True)
    # Remove the history directory that mpet creates.
    try:
        os.rmdir(osp.join(dirDict["suite"], "history"))
    except OSError as exception:
        if exception.errno != errno.ENOENT:
            raise
    return

def compare_with_ref(runInfo, dirDict, tol=1e-7):
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
            if np.max(diffMat) > tol:
                failList.append(testStr)
    return failList

def show_fails(failList):
    for fail in failList:
        print (fail + " differs from the reference outputs!")
    return

def main(compareDir):
    IO = mpetParamsIO.mpetIO()
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
            "test001" : test001,
            "test002" : test002,
            "test003" : test003,
            "test004" : test004,
            "test005" : test005,
            "test006" : test006,
            "test007" : test007,
            "test008" : test008,
            "test009" : test009,
            "test010" : test010,
            "test011" : test011,
            "test012" : test012,
            "test013" : test013,
            "test014" : test014,
            "test015" : test015,
            "test016" : test016,
            "test017" : test017,
            }

    if compareDir is None:
        os.makedirs(dirDict["out"])
        run_test_sims(IO, runInfo, dirDict)
    else:
        dirDict["out"] = compareDir
    failList = compare_with_ref(runInfo, dirDict, tol=1e-7)

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
