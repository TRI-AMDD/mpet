import os.path as osp

import matplotlib.pyplot as plt
import numpy as np

import mpet.plot.plot_data as pd


def vt(testDir, dirDict):
    xyCmp(testDir, dirDict, "vt", "Time", "Voltage",
          "V compare", "vt")


def curr(testDir, dirDict):
    xyCmp(testDir, dirDict, "curr", "Time", "Current",
          "Current compare", "curr")


def elytecons(testDir, dirDict):
    xyCmp(testDir, dirDict, "elytecons", "Time", "Avg elyte conc",
          "Elyte mass cons. -- new min,max:", "elytecons")


def soc(testDir, dirDict, trode):
    xyCmp(testDir, dirDict, "soc_" + trode, "Time", "Filling Frac",
          "Filling fraction " + trode + " compare", "soc_" + trode)


def elytecf(testDir, dirDict):
    xyCmp(testDir, dirDict, "elytecf", "Batt. Posn.", "Elyte conc",
          "Final elyte c profile compare", "elytecf")


def elytepf(testDir, dirDict):
    xyCmp(testDir, dirDict, "elytepf", "Batt. Posn.", "Elyte phi",
          "Final elyte potential profile compare", "elytepf")


def elyteif(testDir, dirDict):
    xyCmp(testDir, dirDict, "elyteif", "Batt. Posn.", "Elyte curr dens",
          "Final elyte curr dens profile compare", "elyteif")


def elytedivif(testDir, dirDict):
    xyCmp(testDir, dirDict, "elytedivif", "Batt. Posn.",
          "div(elyte curr dens)",
          "Final div(elyte curr dens) profile compare", "elytedivif")


def cbarLine(testDir, dirDict, trode):
    xyPartsCmp(testDir, dirDict, "cbar_" + trode, "Time", "cbar",
               "cbar " + trode + " compare", "cbarLine_" + trode)


def bulkpf(testDir, dirDict, trode):
    xyCmp(testDir, dirDict, "bulkpf_" + trode, "Electrode Posn.", "phi",
          "phi of bulk electrode " + trode + " compare", "bulkp_" + trode)


def xyCmp(testDir, dirDict, ptype, xlbl, ylbl, ttl, fname):
    testName = osp.basename(testDir)
    newSimOutDir = osp.join(testDir, "sim_output")
    oldSimOutDir = osp.join(dirDict["refs"], testName, "sim_output")
    plotsDir = dirDict["plots"]
    kwargs = {"print_flag": False, "save_flag": False, "data_only": True}
    new_xy = pd.show_data(newSimOutDir, plot_type=ptype, **kwargs)
    old_xy = pd.show_data(oldSimOutDir, plot_type=ptype, **kwargs)
    scl = 1.3
    fig, ax = plt.subplots(figsize=(scl*6, scl*4))
    ax.plot(new_xy[0], new_xy[1], label="new")
    ax.plot(old_xy[0], old_xy[1], linestyle="--", label="ref")
    ax.legend(loc="best")
    if ptype == "elytecons":
        ttl += "\n{minval:1.6f},{maxval:1.6f}".format(
            minval=np.min(new_xy[1]),
            maxval=np.max(new_xy[1]))
    ax.set_title(ttl)
    ax.set_xlabel(xlbl)
    ax.set_ylabel(ylbl)
    fname = osp.join(plotsDir, testName + "_" + fname)
    fig.savefig(fname, bbox_inches="tight")
    plt.close('all')  # no idea why plt.close(fig) leaves figs open...


def xyPartsCmp(testDir, dirDict, ptype, xlbl, ylbl, ttl, fname):
    testName = osp.basename(testDir)
    newSimOutDir = osp.join(testDir, "sim_output")
    oldSimOutDir = osp.join(dirDict["refs"], testName, "sim_output")
    plotsDir = dirDict["plots"]
    kwargs = {"print_flag": False, "save_flag": False, "data_only": True}
    new_dataDict = pd.show_data(newSimOutDir, plot_type=ptype, **kwargs)
    old_dataDict = pd.show_data(oldSimOutDir, plot_type=ptype, **kwargs)
    new_t = pd.show_data(newSimOutDir, plot_type="vt", **kwargs)[0]
    old_t = pd.show_data(oldSimOutDir, plot_type="vt", **kwargs)[0]
    trode = ptype[-1]
    new_data = new_dataDict[trode]
    old_data = old_dataDict[trode]
    Nt, Nv, Np = new_data.shape
    fig, axs = plt.subplots(Np, Nv, figsize=(6, 4))
    axs = np.reshape(axs, (Np, Nv))
    for i in range(Np):
        for j in range(Nv):
            ax = axs[i, j]
            ax.plot(new_t, new_data[:,j,i], label="new")
            ax.plot(old_t, old_data[:,j,i], linestyle="--", label="ref")
    fig.suptitle(ttl)
    fname = osp.join(plotsDir, testName + "_" + fname)
    fig.savefig(fname, bbox_inches="tight")
    plt.close('all')  # no idea why plt.close(fig) leaves figs open...
