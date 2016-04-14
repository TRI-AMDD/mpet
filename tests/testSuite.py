import os
import sys
import shutil
import time
import ConfigParser

mpetdir = os.path.join(os.environ["HOME"], "docs", "bazantgroup", "mpet")
sys.path.append(mpetdir)
import mpet
import mpetParamsIO

def readConfig(inFile):
    P = ConfigParser.RawConfigParser()
    P.optionxform = str
    P.read(inFile)
    return P

def mod_sys_cfg(IO, psys):
    P = readConfig(psys)
    P.set("Sim Params", "Vmax", "1e10")
    P.set("Sim Params", "Vmin", "-1e10")
    P.set("Sim Params", "tsteps", "200")
    P.set("Sim Params", "relTol", "1e-6")
    P.set("Sim Params", "absTol", "1e-6")
    P.set("Sim Params", "T", "298")
    P.set("Sim Params", "randomSeed", "true")
    P.set("Sim Params", "Rser", "0.")
    P.set("Sim Params", "Nvol_c", "1")
    P.set("Sim Params", "Nvol_s", "0")
    P.set("Sim Params", "Nvol_a", "0")
    P.set("Sim Params", "Npart_c", "1")
    P.set("Sim Params", "Npart_a", "1")
    P.set("Electrodes", "k0_foil", "1e3")
    P.set("Particles", "mean_c", "100e-9")
    P.set("Particles", "stddev_c", "0e-9")
    P.set("Particles", "mean_a", "100e-9")
    P.set("Particles", "stddev_a", "0e-9")
    P.set("Particles", "cs0_c", "0.01")
    P.set("Particles", "cs0_a", "0.99")
    P.set("Conductivity", "simBulkCond_c", "false")
    P.set("Conductivity", "simBulkCond_a", "false")
    P.set("Conductivity", "simPartCond_c", "false")
    P.set("Conductivity", "simPartCond_a", "false")
    P.set("Geometry", "L_c", "50e-6")
    P.set("Geometry", "L_s", "25e-6")
    P.set("Geometry", "L_a", "50e-6")
    P.set("Geometry", "P_L_c", "0.69")
    P.set("Geometry", "P_L_a", "0.95")
    P.set("Geometry", "poros_c", "0.3")
    P.set("Geometry", "poros_s", "0.5")
    P.set("Geometry", "poros_a", "0.3")
    P.set("Electrolytes", "c0", "1000")
    P.set("Electrolytes", "zp", "1")
    P.set("Electrolytes", "zm", "-1")
    P.set("Electrolytes", "nup", "1")
    P.set("Electrolytes", "num", "1")
    P.set("Electrolytes", "elyteModelType", "dilute")
    P.set("Electrolytes", "Dp", "2.2e-10")
    P.set("Electrolytes", "Dm", "2.94e-10")
    IO.writeConfigFile(P, psys)
    return

def test01(IO, psys, ptrodes, run=True):
    # Make a temporary system config file that brings in the right
    # electrode config file
    P_s = readConfig(psys)
    P_s.set("Sim Params", "profileType", "CC")
    P_s.set("Sim Params", "Crate", "1e-2")
    IO.writeConfigFile(P_s, psys)
    P = readConfig(ptrodes["c"])
    P.set("Particles", "type", "ACR")
    P.set("Particles", "discretization", "1e-9")
    P.set("Particles", "shape", "1e-9")
    P.set("Particles", "thickness", "20e-9")
    P.set("Material", "muRfunc", "LiFePO4")
    P.set("Material", "logPad", "false")
    P.set("Material", "noise", "false")
    P.set("Material", "Omega_a", "1.8560e-20")
    P.set("Material", "kappa", "5.0148e-10")
    P.set("Material", "B", "0.1916e9")
    P.set("Material", "rho_s", "1.3793e28")
    P.set("Material", "cwet", "0.98")
    P.set("Reactions", "rxnType", "BV")
    P.set("Reactions", "k0", "1.6e-1")
    P.set("Reactions", "alpha", "0.5")
    IO.writeConfigFile(P, ptrodes["c"])
    P_s, P_e = mpetParamsIO.getConfigs(psys)
    if run:
        mpet.main(psys, keepArchive=False)
    return

def main():
    IO = mpetParamsIO.mpetIO()
    # Get the default configs
    outdir = os.path.join(os.getcwd(),
                          time.strftime("%Y%m%d_%H%M%S", time.localtime()))
    psysName, ptrodeName = "params_system.cfg", "params_electrodes.cfg"
    psys = os.path.join(outdir, psysName)
    ptrodes = {"c" : os.path.join(outdir, "params_c.cfg"),
               "a" : os.path.join(outdir, "params_a.cfg")}
    defaultConfigDir = os.path.join(mpetdir, "configDefaults")
    shutil.copy(os.path.join(defaultConfigDir, psysName), outdir)
    shutil.copy(os.path.join(defaultConfigDir, ptrodeName), ptrodes["c"])
    shutil.copy(os.path.join(defaultConfigDir, ptrodeName), ptrodes["a"])

    # Make an output directory for each test
    test01(psys, ptrodes)
    return

if __name__ == "__main__":
    main()
