import os
import os.path as osp
import sys
import shutil
import time
import ConfigParser

mpetdir = osp.join(os.environ["HOME"], "docs", "bazantgroup", "mpet")
sys.path.append(mpetdir)
import mpet
import mpetParamsIO

def test001(IO, testDir, baseConfigDir, simOutDir, run=True):
    # Make a temporary system config file that brings in the right
    # electrode config file
    shutil.copy(osp.join(baseConfigDir, "params_system.cfg"), testDir)
    shutil.copy(osp.join(baseConfigDir, "params_c.cfg"), testDir)
    psys = osp.join(testDir, "params_system.cfg")
    ptrode = osp.join(testDir, "params_c.cfg")
    P_s = mpetParamsIO.getConfig(psys)
    P_s.set("Sim Params", "profileType", "CC")
    P_s.set("Sim Params", "Crate", "1e-2")
    IO.writeConfigFile(P_s, psys)
    P = mpetParamsIO.getConfig(ptrode)
    P.set("Particles", "type", "ACR")
    P.set("Material", "muRfunc", "LiFePO4")
    IO.writeConfigFile(P, ptrode)
    P_s, P_e = IO.getConfigs(psys)
    if run:
        mpet.main(psys, keepArchive=False)
        shutil.move(simOutDir, osp.join(testDir))
    return

def main():
    IO = mpetParamsIO.mpetIO()
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
#            "test002" : {"fn": test002, "runFlag" : True},
            }

    # Make an output directory for each test
    for testStr in runInfo.keys():
        testDir = osp.join(outdir, testStr)
        os.makedirs(testDir)
        runInfo[testStr]["fn"](
            IO, testDir, baseConfigDir, simOutDir,
            run=runInfo[testStr]["runFlag"])
    # Remove the history directory that mpet creates.
    os.rmdir(osp.join(os.getcwd(), "history"))
    return

if __name__ == "__main__":
    main()
