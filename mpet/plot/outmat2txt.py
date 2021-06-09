"""This can be called to convert the default simulation output (.mat file) to csv files."""
import os

import numpy as np
import h5py

import mpet.plot.plot_data as plot_data
import mpet.utils as utils
from mpet.config import constants

# Strings to be used
RowsStr = "Rows correspond to time points (see generalData.txt).\n"
CCStr = "Columns correspond to cell center positions (see discData.txt)."
FCStr = "Columns correspond to face positions (see discData.txt)."
genDataHdr = ("Time [s], Filling fraction of anode, "
              + "Filling fraction of cathode, "
              + "Voltage [V], Current [C-rate], Current [A/m^2], "
              + "Power [W/m^2]")

zeroStr = "Zero is at the anode current collector."
dccpStr = "discretization cell center positions [m]. "
discCCbattery = ("Battery " + dccpStr + zeroStr)
discCCanode = ("Anode " + dccpStr + zeroStr)
discCCsep = ("Separator " + dccpStr + zeroStr)
discCCcathode = ("Cathode " + dccpStr + zeroStr)
discFC = ("Battery discretization face positions [m]. " + zeroStr)
particleIndxExpl = """
Particle output files are indexed such that Vol 0 is closest to the
anode current collector for both electrodes.  Indexing within volumes
is such that Part 0 is closest to the "carbon backbone" if simPartCond
was set to true for that electrode. Otherwise it is in arbitrary
order.
"""
particleDiscExpl = """
Particle discretization info follows.  Lengths and number of
discretization points are provided for each simulated particle.
Meshes are made as a linear space between zero and particle length
with the given number of points.  Rows correspond to different
simulation volumes with the first row being closest to the anode
current collector and the last closest to the cathode current
collector.  Columns represent individual particles within each
simulated volume in no particular order.
"""

elytecHdr = ("Electrolyte Concentrations [M]\n" + RowsStr + CCStr)
elytepHdr = ("Electrolyte Electric Potential [V]\n" + RowsStr + CCStr)
elyteiHdr = ("Electrolyte Current Density [A/m^2]\n" + RowsStr + FCStr)
elytediviHdr = ("Electrolyte Divergence of Current Density [A/m^3]\n"
                + RowsStr + CCStr)

seeDiscStr = "See discData.txt for particle indexing information."
partStr = "partTrode{l}vol{j}part{i}_"
solffStr = "Solid Filling Fractions"
solTail = ("\n" + RowsStr + CCStr + "\n" + seeDiscStr)
solHdr = (solffStr + solTail)
opStr = ", order parameter "
sol1Hdr = (solffStr + opStr + "1" + solTail)
sol2Hdr = (solffStr + opStr + "2" + solTail)
fnameSolL = "sld{l}Vol{j:03d}Part{i:03d}Conc"
fnameSolR = "Data.txt"
fnameSolBase = (fnameSolL + fnameSolR)
fnameSol1Base = (fnameSolL + "1" + fnameSolR)
fnameSol2Base = (fnameSolL + "2" + fnameSolR)

cbarHdrP1 = ("Average particle filling fractions.\n" + RowsStr)
cbarHdrP2 = """Columns correspond to particles specified by simulated
volume index / particle index within volume.
"""
cbarHdrBase = cbarHdrP1 + cbarHdrP2 + seeDiscStr + "\n"

bulkpHdrP2 = """Columns correspond to this electrode's cell center
positions (see discData.txt).
"""
bulkpHdr = ("Bulk electrode electric potential [V]\n" + RowsStr + bulkpHdrP2)
fnameBulkpBase = "bulkPot{l}Data.txt"


def main(indir, genData=True, discData=True, elyteData=True,
         csldData=True, cbarData=True, bulkpData=True):
    config = plot_data.show_data(
        indir, plot_type="params", print_flag=False, save_flag=False,
        data_only=True)
    trodes = config["trodes"]
    CrateCurr = config["1C_current_density"]  # A/m^2
    psd_len_c = config["psd_len"]["c"]
    Nv_c, Np_c = psd_len_c.shape
    dlm = ","

    def get_trode_str(tr):
        return ("Anode" if tr == "a" else "Cathode")
    if "a" in trodes:
        psd_len_a = config["psd_len"]["a"]
        Nv_a, Np_a = psd_len_a.shape
    tVec, vVec = plot_data.show_data(
        indir, plot_type="vt", print_flag=False, save_flag=False,
        data_only=True)
    ntimes = len(tVec)

    if genData:
        ffVec_c = plot_data.show_data(
            indir, plot_type="soc_c", print_flag=False,
            save_flag=False, data_only=True)[1]
        if "a" in trodes:
            ffVec_a = plot_data.show_data(
                indir, plot_type="soc_a", print_flag=False,
                save_flag=False, data_only=True)[1]
        else:
            ffVec_a = np.ones(len(tVec))
        currVec = plot_data.show_data(
            indir, plot_type="curr", print_flag=False,
            save_flag=False, data_only=True)[1]
        powerVec = plot_data.show_data(
            indir, plot_type="power", print_flag=False,
            save_flag=False, data_only=True)[1]
        genMat = np.zeros((ntimes, 7))
        genMat[:,0] = tVec
        genMat[:,1] = ffVec_a
        genMat[:,2] = ffVec_c
        genMat[:,3] = vVec
        genMat[:,4] = currVec
        genMat[:,5] = currVec * CrateCurr
        genMat[:,6] = powerVec
        np.savetxt(os.path.join(indir, "generalData.txt"),
                   genMat, delimiter=dlm, header=genDataHdr)

    if discData:
        cellCentersVec, facesVec = plot_data.show_data(
            indir, plot_type="discData", print_flag=False,
            save_flag=False, data_only=True)
        with open(os.path.join(indir, "discData.txt"), "w") as fo:
            print(discCCbattery, file=fo)
            print(",".join(map(str, cellCentersVec)), file=fo)
            offset = 0
            if "a" in trodes:
                print(file=fo)
                print(discCCanode, file=fo)
                print(",".join(map(str, cellCentersVec[:Nv_a])), file=fo)
                offset = Nv_a
            print(file=fo)
            print(discCCsep, file=fo)
            print(",".join(map(str, cellCentersVec[offset:-Nv_c])), file=fo)
            print(file=fo)
            print(discCCcathode, file=fo)
            print(",".join(map(str, cellCentersVec[-Nv_c:])), file=fo)
            print(file=fo)
            print(discFC, file=fo)
            print(",".join(map(str, facesVec)), file=fo)
            print(file=fo)
            print(particleIndxExpl, file=fo)
            print(particleDiscExpl, file=fo)
            for tr in trodes:
                print(file=fo)
                Trode = get_trode_str(tr)
                print((Trode + " particle sizes [m]"), file=fo)
                for vind in range(config["Nvol"][tr]):
                    print(",".join(map(str, config["psd_len"][tr][vind,:])), file=fo)
                print("\n" + Trode + " particle number of discr. points", file=fo)
                for vind in range(config["Nvol"][tr]):
                    print(",".join(map(str, config["psd_num"][tr][vind,:])), file=fo)

    if elyteData:
        valid_current = True
        # If there wasn't an electrolyte, then a ghost point variable wouldn't have been created,
        # so we'll get a KeyError in attempting to "plot" the electrolyte current density.
        try:
            plot_data.show_data(
                indir, plot_type="elytei", print_flag=False, save_flag=False, data_only=True)
        except KeyError:
            valid_current = False
        elytecMat = plot_data.show_data(
            indir, plot_type="elytec", print_flag=False,
            save_flag=False, data_only=True)[1]
        elytepMat = plot_data.show_data(
            indir, plot_type="elytep", print_flag=False,
            save_flag=False, data_only=True)[1]
        np.savetxt(os.path.join(indir, "elyteConcData.txt"),
                   elytecMat, delimiter=dlm, header=elytecHdr)
        np.savetxt(os.path.join(indir, "elytePotData.txt"),
                   elytepMat, delimiter=dlm, header=elytepHdr)
        if valid_current:
            elyteiMat = plot_data.show_data(
                indir, plot_type="elytei", print_flag=False,
                save_flag=False, data_only=True)[1]
            elytediviMat = plot_data.show_data(
                indir, plot_type="elytedivi", print_flag=False,
                save_flag=False, data_only=True)[1]
            np.savetxt(os.path.join(indir, "elyteCurrDensData.txt"),
                       elyteiMat, delimiter=dlm, header=elyteiHdr)
            np.savetxt(os.path.join(indir, "elyteDivCurrDensData.txt"),
                       elytediviMat, delimiter=dlm, header=elytediviHdr)

    if csldData:
        dataFileName = "output_data"
        dataFile = os.path.join(indir, dataFileName)
        data = utils.open_data_file(dataFile)
        for tr in trodes:
            Trode = get_trode_str(tr)
            type2c = False
            if config[tr, "type"] in constants.one_var_types:
                str_base = partStr + "c"
            elif config[tr, "type"] in constants.two_var_types:
                type2c = True
                str1_base = partStr + "c1"
                str2_base = partStr + "c2"
            for i in range(config["Npart"][tr]):
                for j in range(config["Nvol"][tr]):
                    if type2c:
                        sol1 = str1_base.format(l=tr, i=i, j=j)
                        sol2 = str2_base.format(l=tr, i=i, j=j)
                        datay1 = utils.get_dict_key(data, sol1)
                        datay2 = utils.get_dict_key(data, sol2)
                        filename1 = fnameSol1Base.format(l=Trode, i=i, j=j)
                        filename2 = fnameSol2Base.format(l=Trode, i=i, j=j)
                        np.savetxt(os.path.join(indir, filename1), datay1,
                                   delimiter=dlm, header=sol1Hdr)
                        np.savetxt(os.path.join(indir, filename2), datay2,
                                   delimiter=dlm, header=sol2Hdr)
                    else:
                        sol = str_base.format(l=tr, i=i, j=j)
                        datay = utils.get_dict_key(data, sol)
                        filename = fnameSolBase.format(l=Trode, i=i, j=j)
                        np.savetxt(os.path.join(indir, filename), datay,
                                   delimiter=dlm, header=solHdr)

        # close file if it is a h5py file
        if isinstance(data, h5py._hl.files.File):
            data.close()

    if cbarData:
        cbarDict = plot_data.show_data(
            indir, plot_type="cbar_full", print_flag=False,
            save_flag=False, data_only=True)
        for tr in trodes:
            Trode = get_trode_str(tr)
            fname = "cbar{l}Data.txt".format(l=Trode)
            Nv, Np = config["Nvol"][tr], config["Npart"][tr]
            NpartTot = Nv*Np
            cbarMat = np.zeros((ntimes, NpartTot))
            cbarHdr = cbarHdrBase
            partInd = 0
            for i in range(Nv):
                for j in range(Np):
                    cbarHdr += "{i}/{j},".format(j=j, i=i)
                    cbarMat[:,partInd] = cbarDict[tr][:,i,j]
                    partInd += 1
            np.savetxt(os.path.join(indir, fname), cbarMat,
                       delimiter=dlm, header=cbarHdr.rstrip(','))

    if bulkpData:
        if "a" in trodes:
            bulkp_aData = plot_data.show_data(
                indir, plot_type="bulkp_a", print_flag=False,
                save_flag=False, data_only=True)[1]
            fname = fnameBulkpBase.format(l="Anode")
            np.savetxt(os.path.join(indir, fname), bulkp_aData,
                       delimiter=dlm, header=bulkpHdr)
        bulkp_cData = plot_data.show_data(
            indir, plot_type="bulkp_c", print_flag=False,
            save_flag=False, data_only=True)[1]
        fname = fnameBulkpBase.format(l="Cathode")
        np.savetxt(os.path.join(indir, fname), bulkp_cData,
                   delimiter=dlm, header=bulkpHdr)

    return
