import os
import os.path as osp
import sys
import errno
import shutil
import time
import ConfigParser

mpetdir = osp.join(os.environ["HOME"], "docs", "bazantgroup", "mpet")
sys.path.append(mpetdir)
import mpet
import mpetParamsIO as IO

def test001(testDir, baseConfigDir, simOutDir, run=True):
    """ LFP ACR C3 """
    shutil.copy(osp.join(baseConfigDir, "params_system.cfg"), testDir)
    shutil.copy(osp.join(baseConfigDir, "params_c.cfg"), testDir)
    psys = osp.join(testDir, "params_system.cfg")
    ptrode = osp.join(testDir, "params_c.cfg")
    P_s = IO.getConfig(psys)
    P_s.set("Sim Params", "Crate", "1e-2")
    IO.writeConfigFile(P_s, psys)
    P = IO.getConfig(ptrode)
    P.set("Particles", "type", "ACR")
    P.set("Material", "muRfunc", "LiFePO4")
    IO.writeConfigFile(P, ptrode)
    if run:
        mpet.main(psys, keepArchive=False)
        shutil.move(simOutDir, osp.join(testDir))
    return

def test002(testDir, baseConfigDir, simOutDir, run=True):
    """ LFP CHR cylinder """
    shutil.copy(osp.join(baseConfigDir, "params_system.cfg"), testDir)
    shutil.copy(osp.join(baseConfigDir, "params_c.cfg"), testDir)
    psys = osp.join(testDir, "params_system.cfg")
    ptrode = osp.join(testDir, "params_c.cfg")
    P_s = IO.getConfig(psys)
    P_s.set("Sim Params", "Crate", "1e-2")
    IO.writeConfigFile(P_s, psys)
    P = IO.getConfig(ptrode)
    P.set("Particles", "type", "CHR")
    P.set("Particles", "shape", "cylinder")
    P.set("Material", "muRfunc", "LiFePO4")
    P.set("Material", "dgammadc", "5e-30")
    IO.writeConfigFile(P, ptrode)
    if run:
        mpet.main(psys, keepArchive=False)
        shutil.move(simOutDir, osp.join(testDir))
    return

def test003(testDir, baseConfigDir, simOutDir, run=True):
    """ LFP CHR sphere """
    shutil.copy(osp.join(baseConfigDir, "params_system.cfg"), testDir)
    shutil.copy(osp.join(baseConfigDir, "params_c.cfg"), testDir)
    psys = osp.join(testDir, "params_system.cfg")
    ptrode = osp.join(testDir, "params_c.cfg")
    P_s = IO.getConfig(psys)
    P_s.set("Sim Params", "Crate", "1e-2")
    IO.writeConfigFile(P_s, psys)
    P = IO.getConfig(ptrode)
    P.set("Particles", "type", "CHR")
    P.set("Particles", "shape", "sphere")
    P.set("Material", "muRfunc", "LiFePO4")
    P.set("Material", "dgammadc", "-2e-30")
    IO.writeConfigFile(P, ptrode)
    if run:
        mpet.main(psys, keepArchive=False)
        shutil.move(simOutDir, osp.join(testDir))
    return

def test004(testDir, baseConfigDir, simOutDir, run=True):
    """ LFP CHR sphere with noise  """
    shutil.copy(osp.join(baseConfigDir, "params_system.cfg"), testDir)
    shutil.copy(osp.join(baseConfigDir, "params_c.cfg"), testDir)
    psys = osp.join(testDir, "params_system.cfg")
    ptrode = osp.join(testDir, "params_c.cfg")
    P_s = IO.getConfig(psys)
    P_s.set("Sim Params", "Crate", "1e-2")
    IO.writeConfigFile(P_s, psys)
    P = IO.getConfig(ptrode)
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

def test005(testDir, baseConfigDir, simOutDir, run=True):
    """ LFP homog """
    shutil.copy(osp.join(baseConfigDir, "params_system.cfg"), testDir)
    shutil.copy(osp.join(baseConfigDir, "params_c.cfg"), testDir)
    psys = osp.join(testDir, "params_system.cfg")
    ptrode = osp.join(testDir, "params_c.cfg")
    P_s = IO.getConfig(psys)
    P_s.set("Sim Params", "Crate", "1e-2")
    P_s.set("Sim Params", "Nvol_c", "2")
    P_s.set("Sim Params", "Nvol_s", "2")
    IO.writeConfigFile(P_s, psys)
    P = IO.getConfig(ptrode)
    P.set("Particles", "type", "homog")
    P.set("Material", "muRfunc", "LiFePO4")
    IO.writeConfigFile(P, ptrode)
    if run:
        mpet.main(psys, keepArchive=False)
        shutil.move(simOutDir, osp.join(testDir))
    return

def test006(testDir, baseConfigDir, simOutDir, run=True):
    """ LFP homog with logPad, Vmin """
    shutil.copy(osp.join(baseConfigDir, "params_system.cfg"), testDir)
    shutil.copy(osp.join(baseConfigDir, "params_c.cfg"), testDir)
    psys = osp.join(testDir, "params_system.cfg")
    ptrode = osp.join(testDir, "params_c.cfg")
    P_s = IO.getConfig(psys)
    P_s.set("Sim Params", "Crate", "1e-2")
    P_s.set("Sim Params", "Vmin", "2.8e-0")
    P_s.set("Particles", "cs0_c", "0.5")
    IO.writeConfigFile(P_s, psys)
    P = IO.getConfig(ptrode)
    P.set("Particles", "type", "homog")
    P.set("Material", "muRfunc", "LiFePO4")
    IO.writeConfigFile(P, ptrode)
    if run:
        mpet.main(psys, keepArchive=False)
        shutil.move(simOutDir, osp.join(testDir))
    return

def test007(testDir, baseConfigDir, simOutDir, run=True):
    """ LFP homog with logPad, Vmax """
    shutil.copy(osp.join(baseConfigDir, "params_system.cfg"), testDir)
    shutil.copy(osp.join(baseConfigDir, "params_c.cfg"), testDir)
    psys = osp.join(testDir, "params_system.cfg")
    ptrode = osp.join(testDir, "params_c.cfg")
    P_s = IO.getConfig(psys)
    P_s.set("Sim Params", "Crate", "-1e-2")
    P_s.set("Sim Params", "Vmax", "4.0e-0")
    P_s.set("Particles", "cs0_c", "0.5")
    IO.writeConfigFile(P_s, psys)
    P = IO.getConfig(ptrode)
    P.set("Particles", "type", "homog")
    P.set("Material", "muRfunc", "LiFePO4")
    IO.writeConfigFile(P, ptrode)
    if run:
        mpet.main(psys, keepArchive=False)
        shutil.move(simOutDir, osp.join(testDir))
    return

def test008(testDir, baseConfigDir, simOutDir, run=True):
    """ LFP homog_sdn """
    shutil.copy(osp.join(baseConfigDir, "params_system.cfg"), testDir)
    shutil.copy(osp.join(baseConfigDir, "params_c.cfg"), testDir)
    psys = osp.join(testDir, "params_system.cfg")
    ptrode = osp.join(testDir, "params_c.cfg")
    P_s = IO.getConfig(psys)
    P_s.set("Sim Params", "Crate", "1e-2")
    P_s.set("Sim Params", "Nvol_c", "3")
    P_s.set("Sim Params", "Npart_c", "3")
    P_s.set("Particles", "stddev_c", "25e-9")
    P_s.set("Sim Params", "capFrac", "0.6")
    IO.writeConfigFile(P_s, psys)
    P = IO.getConfig(ptrode)
    P.set("Particles", "type", "homog_sdn")
    P.set("Material", "muRfunc", "LiFePO4")
    IO.writeConfigFile(P, ptrode)
    if run:
        mpet.main(psys, keepArchive=False)
        shutil.move(simOutDir, osp.join(testDir))
    return

def test009(testDir, baseConfigDir, simOutDir, run=True):
    """ Graphite-2param homog """
    shutil.copy(osp.join(baseConfigDir, "params_system.cfg"), testDir)
    shutil.copy(osp.join(baseConfigDir, "params_c.cfg"), testDir)
    psys = osp.join(testDir, "params_system.cfg")
    ptrode = osp.join(testDir, "params_c.cfg")
    P_s = IO.getConfig(psys)
    P_s.set("Sim Params", "Crate", "1e-2")
    IO.writeConfigFile(P_s, psys)
    P = IO.getConfig(ptrode)
    P.set("Particles", "type", "homog2")
    P.set("Material", "muRfunc", "LiC6")
    IO.writeConfigFile(P, ptrode)
    if run:
        mpet.main(psys, keepArchive=False)
        shutil.move(simOutDir, osp.join(testDir))
    return

def test010(testDir, baseConfigDir, simOutDir, run=True):
    """ Graphite-2param CHR cylinder """
    shutil.copy(osp.join(baseConfigDir, "params_system.cfg"), testDir)
    shutil.copy(osp.join(baseConfigDir, "params_c.cfg"), testDir)
    psys = osp.join(testDir, "params_system.cfg")
    ptrode = osp.join(testDir, "params_c.cfg")
    P_s = IO.getConfig(psys)
    P_s.set("Sim Params", "Crate", "1e-2")
    IO.writeConfigFile(P_s, psys)
    P = IO.getConfig(ptrode)
    P.set("Particles", "type", "CHR2")
    P.set("Particles", "shape", "cylinder")
    P.set("Material", "muRfunc", "LiC6")
    IO.writeConfigFile(P, ptrode)
    if run:
        mpet.main(psys, keepArchive=False)
        shutil.move(simOutDir, osp.join(testDir))
    return

def test011(testDir, baseConfigDir, simOutDir, run=True):
    """ Graphite-2param CHR sphere """
    shutil.copy(osp.join(baseConfigDir, "params_system.cfg"), testDir)
    shutil.copy(osp.join(baseConfigDir, "params_c.cfg"), testDir)
    psys = osp.join(testDir, "params_system.cfg")
    ptrode = osp.join(testDir, "params_c.cfg")
    P_s = IO.getConfig(psys)
    P_s.set("Sim Params", "Crate", "1e-2")
    P_s.set("Sim Params", "relTol", "1e-7")
    P_s.set("Sim Params", "absTol", "1e-7")
    IO.writeConfigFile(P_s, psys)
    P = IO.getConfig(ptrode)
    P.set("Particles", "type", "CHR2")
    P.set("Particles", "shape", "sphere")
    P.set("Material", "muRfunc", "LiC6")
    IO.writeConfigFile(P, ptrode)
    if run:
        mpet.main(psys, keepArchive=False)
        shutil.move(simOutDir, osp.join(testDir))
    return

def test012(testDir, baseConfigDir, simOutDir, run=True):
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
    P_s = IO.getConfig(psys)
    P_s.set("Sim Params", "Crate", "1e-2")
    P_s.set("Sim Params", "capFrac", "0.67")
    P_s.set("Sim Params", "Nvol_c", "2")
    P_s.set("Sim Params", "Nvol_a", "2")
    P_s.set("Sim Params", "Vmin", "2e0")
    P_s.set("Particles", "cs0_c", "0.2")
    P_s.set("Particles", "cs0_a", "0.495")
    P_s.set("Electrolyte", "elyteModelType", "SM")
    IO.writeConfigFile(P_s, psys)
    P = IO.getConfig(ptrodec)
    P.set("Particles", "type", "homog")
    P.set("Particles", "shape", "cylinder")
    P.set("Material", "muRfunc", "LiMn2O4_ss2")
    P.set("Reactions", "rxnType", "BV_mod01")
    IO.writeConfigFile(P, ptrodec)
    P = IO.getConfig(ptrodea)
    P.set("Particles", "shape", "sphere")
    P.set("Material", "muRfunc", "LiC6_coke_ss2")
    P.set("Reactions", "rxnType", "BV_mod02")
    IO.writeConfigFile(P, ptrodea)
    if run:
        mpet.main(psys, keepArchive=False)
        shutil.move(simOutDir, osp.join(testDir))
    return

def test013(testDir, baseConfigDir, simOutDir, run=True):
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
    P_s = IO.getConfig(psys)
    P_s.set("Sim Params", "profileType", "CV")
    P_s.set("Sim Params", "Vset", "3.8")
    P_s.set("Sim Params", "Nvol_c", "2")
    P_s.set("Sim Params", "Nvol_s", "2")
    P_s.set("Sim Params", "Nvol_a", "2")
    P_s.set("Particles", "cs0_c", "0.2")
    P_s.set("Particles", "cs0_a", "0.95")
    P_s.set("Electrolyte", "elyteModelType", "dilute")
    IO.writeConfigFile(P_s, psys)
    P = IO.getConfig(ptrodec)
    P.set("Particles", "type", "homog")
    P.set("Particles", "shape", "sphere")
    P.set("Material", "muRfunc", "LiMn2O4_ss2")
    P.set("Reactions", "rxnType", "Marcus")
    IO.writeConfigFile(P, ptrodec)
    P = IO.getConfig(ptrodea)
    P.set("Particles", "shape", "cylinder")
    P.set("Material", "muRfunc", "testIS_ss")
    P.set("Reactions", "rxnType", "BV_raw")
    IO.writeConfigFile(P, ptrodea)
    if run:
        mpet.main(psys, keepArchive=False)
        shutil.move(simOutDir, osp.join(testDir))
    return

def test014(testDir, baseConfigDir, simOutDir, run=True):
    """ LFP homog with CCsegments, MHC, Rser """
    shutil.copy(osp.join(baseConfigDir, "params_system.cfg"), testDir)
    shutil.copy(osp.join(baseConfigDir, "params_c.cfg"), testDir)
    psys = osp.join(testDir, "params_system.cfg")
    ptrode = osp.join(testDir, "params_c.cfg")
    P_s = IO.getConfig(psys)
    P_s.set("Sim Params", "profileType", "CCsegments")
    P_s.set("Sim Params", "segments",
            "[(1., 25), (-2., 10), (0., 30)]")
    P_s.set("Sim Params", "Rser", "1e-3")
    IO.writeConfigFile(P_s, psys)
    P = IO.getConfig(ptrode)
    P.set("Particles", "type", "homog")
    P.set("Material", "muRfunc", "LiFePO4")
    P.set("Reactions", "rxnType", "MHC")
    IO.writeConfigFile(P, ptrode)
    if run:
        mpet.main(psys, keepArchive=False)
        shutil.move(simOutDir, osp.join(testDir))
    return

def test015(testDir, baseConfigDir, simOutDir, run=True):
    """ testRS homog with CVsegments, bulkCond, partCond """
    shutil.copy(osp.join(baseConfigDir, "params_system.cfg"), testDir)
    shutil.copy(osp.join(baseConfigDir, "params_c.cfg"), testDir)
    psys = osp.join(testDir, "params_system.cfg")
    ptrode = osp.join(testDir, "params_c.cfg")
    P_s = IO.getConfig(psys)
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
    P = IO.getConfig(ptrode)
    P.set("Particles", "type", "homog")
    P.set("Material", "muRfunc", "testRS")
    P.set("Reactions", "rxnType", "BV_raw")
    IO.writeConfigFile(P, ptrode)
    if run:
        mpet.main(psys, keepArchive=False)
        shutil.move(simOutDir, osp.join(testDir))
    return

def test016(testDir, baseConfigDir, simOutDir, run=True):
    """ test CC continuation """
    shutil.copy(osp.join(baseConfigDir, "params_system.cfg"), testDir)
    shutil.copy(osp.join(baseConfigDir, "params_c.cfg"), testDir)
    psys = osp.join(testDir, "params_system.cfg")
    ptrode = osp.join(testDir, "params_c.cfg")
    P_s = IO.getConfig(psys)
    P_s.set("Sim Params", "Crate", "1e-2")
    P_s.set("Sim Params", "Nvol_c", "3")
    P_s.set("Sim Params", "Npart_c", "3")
    P_s.set("Sim Params", "capFrac", "0.34")
    test008dir = str(osp.join(testDir, "..", "test008", "sim_output"))
    P_s.set("Sim Params", "prevDir", test008dir)
    IO.writeConfigFile(P_s, psys)
    P = IO.getConfig(ptrode)
    P.set("Particles", "type", "homog_sdn")
    P.set("Material", "muRfunc", "LiFePO4")
    IO.writeConfigFile(P, ptrode)
    if run:
        mpet.main(psys, keepArchive=False)
        shutil.move(simOutDir, osp.join(testDir))
    return

def test017(testDir, baseConfigDir, simOutDir, run=True):
    """ test CV continuation """
    shutil.copy(osp.join(baseConfigDir, "params_system.cfg"), testDir)
    shutil.copy(osp.join(baseConfigDir, "params_c.cfg"), testDir)
    psys = osp.join(testDir, "params_system.cfg")
    ptrode = osp.join(testDir, "params_c.cfg")
    P_s = IO.getConfig(psys)
    P_s.set("Sim Params", "profileType", "CV")
    P_s.set("Sim Params", "Vset", "3.45")
    P_s.set("Sim Params", "tend", "3e3")
    P_s.set("Sim Params", "Nvol_c", "3")
    P_s.set("Sim Params", "Npart_c", "3")
    test008dir = str(osp.join(testDir, "..", "test008", "sim_output"))
    P_s.set("Sim Params", "prevDir", test008dir)
    IO.writeConfigFile(P_s, psys)
    P = IO.getConfig(ptrode)
    P.set("Particles", "type", "homog_sdn")
    P.set("Material", "muRfunc", "LiFePO4")
    IO.writeConfigFile(P, ptrode)
    if run:
        mpet.main(psys, keepArchive=False)
        shutil.move(simOutDir, osp.join(testDir))
    return

def test018(testDir, baseConfigDir, simOutDir, run=True):
    """ Like test014, LFP homog with CCsegments, BV, Rfilm, Rfilm_foil """
    shutil.copy(osp.join(baseConfigDir, "params_system.cfg"), testDir)
    shutil.copy(osp.join(baseConfigDir, "params_c.cfg"), testDir)
    psys = osp.join(testDir, "params_system.cfg")
    ptrode = osp.join(testDir, "params_c.cfg")
    P_s = IO.getConfig(psys)
    P_s.set("Sim Params", "profileType", "CCsegments")
    P_s.set("Sim Params", "segments",
            "[(1., 25), (-2., 10), (0., 30)]")
    P_s.set("Electrodes", "k0_foil", "3e+0")
    P_s.set("Electrodes", "Rfilm_foil", "3e+0")
    IO.writeConfigFile(P_s, psys)
    P = IO.getConfig(ptrode)
    P.set("Particles", "type", "homog")
    P.set("Material", "muRfunc", "LiFePO4")
    P.set("Reactions", "Rfilm", "1e+1")
    IO.writeConfigFile(P, ptrode)
    if run:
        mpet.main(psys, keepArchive=False)
        shutil.move(simOutDir, osp.join(testDir))
    return

def main():
    # Get the default configs
    suiteDir = osp.dirname(osp.abspath(__file__))
    simOutDir = osp.join(os.getcwd(), "sim_output")
    outdir = osp.join(suiteDir,
                      time.strftime("%Y%m%d_%H%M%S", time.localtime()))
    os.makedirs(outdir)
    baseConfigDir = osp.join(suiteDir, "baseConfigs")

    # Dictionary containing info about the tests to run
    # Identifier strings are associated with functions to call and
    # whether to run that particular test.
    runInfo = {
            "test001" : {"fn": test001, "runFlag" : True},
            "test002" : {"fn": test002, "runFlag" : True},
            "test003" : {"fn": test003, "runFlag" : True},
            "test004" : {"fn": test004, "runFlag" : True},
            "test005" : {"fn": test005, "runFlag" : True},
            "test006" : {"fn": test006, "runFlag" : True},
            "test007" : {"fn": test007, "runFlag" : True},
            "test008" : {"fn": test008, "runFlag" : True},
            "test009" : {"fn": test009, "runFlag" : True},
            "test010" : {"fn": test010, "runFlag" : True},
            "test011" : {"fn": test011, "runFlag" : True},
            "test012" : {"fn": test012, "runFlag" : True},
            "test013" : {"fn": test013, "runFlag" : True},
            "test014" : {"fn": test014, "runFlag" : True},
            "test015" : {"fn": test015, "runFlag" : True},
            "test016" : {"fn": test016, "runFlag" : True},
            "test017" : {"fn": test017, "runFlag" : True},
            "test018" : {"fn": test018, "runFlag" : True},
            }

    # Make an output directory for each test
    for testStr in sorted(runInfo.keys()):
        testDir = osp.join(outdir, testStr)
        os.makedirs(testDir)
        runInfo[testStr]["fn"](
            testDir, baseConfigDir, simOutDir,
            run=runInfo[testStr]["runFlag"])
    # Remove the history directory that mpet creates.
    try:
        os.rmdir(osp.join(suiteDir, "history"))
    except OSError as exception:
        if exception.errno != errno.ENOENT:
            raise
    return

if __name__ == "__main__":
    main()
