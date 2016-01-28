import sys
import os

import numpy as np
import scipy.io as sio

import plot_data

def main(indir, genData=True, discData=True, elyteData=True,
        csldData=True, cbarData=True, bulkpData=True):
    # First file:
    # time, ff_a, ff_c, voltage, current
    # Second file:
    # discretization along battery
    # discretization information for particles
    # Third file:
    # electrolyte profiles
    # Fourth file:
    # solid particle profiles?
    ndD_s, dD_s, ndD_e, dD_e = plot_data.show_data(indir,
            plot_type="params", print_flag=False,
            save_flag=False, data_only=True)
    limtrode = ("c" if ndD_s["z"] < 1 else "a")
    trodes = ndD_s["trodes"]
    CrateCurr = dD_e[limtrode]["cap"] / 3600. # A/m^2
    psd_len_c = dD_s["psd_len"]["c"]
    psd_num_c = ndD_s["psd_num"]["c"]
    Nv_c, Np_c = psd_len_c.shape
    if "a" in trodes:
        psd_len_a = dD_s["psd_len"]["a"]
        psd_num_a = ndD_s["psd_num"]["a"]
        Nv_a, Np_a = psd_len_a.shape
    RowsStr = "Rows correspond to time points (see generalData.txt).\n"
    CCStr = "Columns correspond to cell center positions (see discData.txt)."
    FCStr = "Columns correspond to face positions (see discData.txt)."
    tVec, vVec = plot_data.show_data(indir, plot_type="vt",
            print_flag=False, save_flag=False, data_only=True)
    ntimes = len(tVec)

    if genData:
        ffVec_c = plot_data.show_data(indir, plot_type="soc_c",
                print_flag=False, save_flag=False, data_only=True)[1]
        if "a" in trodes:
            ffVec_a = plot_data.show_data(indir, plot_type="soc_a",
                    print_flag=False, save_flag=False,
                    data_only=True)[1]
        else:
            ffVec_a = np.ones(len(tVec))
        currVec = plot_data.show_data(indir, plot_type="curr",
                print_flag=False, save_flag=False, data_only=True)[1]
        genMat = np.zeros((ntimes, 6))
        genMat[:, 0] = tVec
        genMat[:, 1] = ffVec_a
        genMat[:, 2] = ffVec_c
        genMat[:, 3] = vVec
        genMat[:, 4] = currVec
        genMat[:, 5] = currVec * CrateCurr
        genStr = ("Time [s], Filling fraction of anode, " +
                  "Filling fraction of cathode, " +
                  "Voltage [V], Current [C-rate], " +
                  "Current [A/m^2]")
        np.savetxt(os.path.join(indir, "generalData.txt"),
                genMat, delimiter=",", header=genStr)

    if discData:
        cellCentersVec, facesVec = plot_data.show_data(indir,
                plot_type="discData",
                print_flag=False, save_flag=False, data_only=True)
        with open(os.path.join(indir, "discData.txt"), "w") as fo:
            print >> fo, "Battery discretization cell center positions [m]"
            print >> fo, "Zero is at the anode current collector."
            print >> fo, ",".join(map(str, cellCentersVec))
            offset = 0
            if "a" in trodes:
                print >> fo
                print >> fo, "Anode discretization cell center positions [m]"
                print >> fo, ",".join(map(str, cellCentersVec[:Nv_a]))
                offset = Nv_a
            print >> fo
            print >> fo, "Separator discretization cell center positions [m]"
            print >> fo, ",".join(map(str, cellCentersVec[offset:-Nv_c]))
            print >> fo
            print >> fo, "Cathode discretization cell center positions [m]"
            print >> fo, ",".join(map(str, cellCentersVec[-Nv_c:]))
            print >> fo
            print >> fo, "Battery discretization face positions [m]"
            print >> fo, ",".join(map(str, facesVec))
            print >> fo
            print >> fo, ("Particle discretization info follows. " +
                    "Lengths and number of discretization points " +
                    "are provided for each simulated particle.\n" +
                    "Meshes are made as a linear space between " +
                    "zero and particle length with the given " +
                    "number of points.\n" +
                    "Rows correspond to different simulation volumes " +
                    "with the first row being closest to the " +
                    "anode current collector and the last " +
                    "closest to the cathode current collector.\n"
                    "Columns represent individual particles within " +
                    "each simulated volume in no particular order.")
            print >> fo
            if "a" in trodes:
                print >> fo, ("Anode particle sizes [m]")
                for vind in range(Nv_a):
                    print >> fo, ",".join(map(str,
                        psd_len_a[vind, :]))
                print >> fo, ("Anode particle number of discr. points")
                for vind in range(Nv_a):
                    print >> fo, ",".join(map(str,
                        psd_num_a[vind, :]))

            print >> fo
            print >> fo, ("Cathode particle sizes [m]")
            for vind in range(Nv_c):
                print >> fo, ",".join(map(str,
                    psd_len_c[vind, :]))
            print >> fo, ("Cathode particle number of discr. points")
            for vind in range(Nv_c):
                print >> fo, ",".join(map(str,
                    psd_num_c[vind, :]))

    if elyteData:
        elytecMat = plot_data.show_data(indir,
                plot_type="elytec",
                print_flag=False, save_flag=False, data_only=True)[1]
        elytepMat = plot_data.show_data(indir,
                plot_type="elytep",
                print_flag=False, save_flag=False, data_only=True)[1]
        elyteiMat = plot_data.show_data(indir,
                plot_type="elytei",
                print_flag=False, save_flag=False, data_only=True)[1]
        elytediviMat = plot_data.show_data(indir,
                plot_type="elytedivi",
                print_flag=False, save_flag=False, data_only=True)[1]
        elytecStr = ("Electrolyte Concentrations [M]\n" +
                RowsStr + CCStr)
        elytepStr = ("Electrolyte Electric Potential [V]\n" +
                RowsStr + CCStr)
        elyteiStr = ("Electrolyte Current Density [A/m^2]\n" +
                RowsStr + FCStr)
        elytediviStr = ("Electrolyte Divergence of Current " +
                "Density [A/m^3]\n" +
                RowsStr + CCStr)
        np.savetxt(os.path.join(indir, "elyteConcData.txt"),
                elytecMat, delimiter=',', header=elytecStr)
        np.savetxt(os.path.join(indir, "elytePotData.txt"),
                elytepMat, delimiter=',', header=elytepStr)
        np.savetxt(os.path.join(indir, "elyteCurrDensData.txt"),
                elyteiMat, delimiter=',', header=elyteiStr)
        np.savetxt(os.path.join(indir, "elyteDivCurrDensData.txt"),
                elytediviMat, delimiter=',', header=elytediviStr)

    if csldData:
        ttl_fmt = "% = {perc:2.1f}"
        dataFileName = "output_data.mat"
        dataFile = os.path.join(indir, dataFileName)
        data = sio.loadmat(dataFile)
        partStr = "partTrode{l}vol{j}part{i}."
        type2c = False
        for l in trodes:
            Trode = ("Anode" if l == "a" else "Cathode")
            if ndD_e[l]["type"] in ndD_s["1varTypes"]:
                str_base = partStr + "c"
            elif ndD_e[l]["type"] in ndD_s["2varTypes"]:
                type2c = True
                str1_base = partStr + "c1"
                str2_base = partStr + "c2"
            for i in range(ndD_s["Npart"][l]):
                for j in range(ndD_s["Nvol"][l]):
                    if type2c:
                        sol1 = str1_base.format(l=l, i=i, j=j)
                        sol2 = str2_base.format(l=l, i=i, j=j)
                        datay1 = data[sol1]
                        datay2 = data[sol2]
                        sol1Str = ("Solid Filling Fractions, order parameter 1\n" +
                                RowsStr + CCStr)
                        sol2Str = ("Solid Filling Fractions, order parameter 2\n" +
                                RowsStr + CCStr)
                        filename1 = "sld{l}Vol{j:03d}Part{i:03d}Conc1Data.txt".format(
                                l=Trode,i=i,j=j)
                        filename2 = "sld{l}Vol{j:03d}Part{i:03d}Conc2Data.txt".format(
                                l=Trode,i=i,j=j)
                        np.savetxt(os.path.join(indir, filename1), datay1,
                                delimiter=",", header=sol1Str)
                        np.savetxt(os.path.join(indir, filename2), datay2,
                                delimiter=",", header=sol2Str)
                    else:
                        sol = str_base.format(l=l, i=i, j=j)
                        datay = data[sol]
                        solStr = ("Solid Filling Fractions\n" +
                                RowsStr + CCStr)
                        filename = "sld{l}Vol{j:03d}Part{i:03d}ConcData.txt".format(
                                l=Trode,i=i,j=j)
                        np.savetxt(os.path.join(indir, filename), datay,
                                delimiter=",", header=solStr)

    if cbarData:
        cbarDict = plot_data.show_data(indir, plot_type="cbar_full",
                print_flag=False, save_flag=False, data_only=True)
        for l in trodes:
            Trode = ("Anode" if l == "a" else "Cathode")
            fname = "cbar{l}Data.txt".format(l=Trode)
            Nv, Np = ndD_s["Nvol"][l], ndD_s["Npart"][l]
            NpartTot = Nv*Np
            cbarMat = np.zeros((ntimes, NpartTot))
            cbarStr = ("Average particle filling fractions.\n" +
                    RowsStr +
                    "Columns correspond to particles specified by\n" +
                    "Simulated volume index / particle index " +
                    "within volume. " +
                    "See discData.txt for more details.\n")
            partInd = 0
            for i in range(Nv):
                for j in range(Np):
                    cbarStr += "{i}/{j},".format(j=j,i=i)
                    cbarMat[:, partInd] = cbarDict[l][:, i, j]
                    partInd += 1
            np.savetxt(os.path.join(indir, fname), cbarMat,
                    delimiter=",", header=cbarStr)

    if bulkpData:
        CCTrodeStr = ("Columns correspond to this electrode's cell " +
                "center positions (see discData.txt)")
        header = ("Bulk electrode electric potential [V]\n"
                + RowsStr + CCTrodeStr)
        if "a" in trodes:
            bulkp_aData = plot_data.show_data(indir, plot_type="bulkp_a",
                print_flag=False, save_flag=False, data_only=True)[1]
            fname = "bulkPotAnodeData.txt"
            np.savetxt(os.path.join(indir, fname), bulkp_aData,
                    delimiter=",", header=header)
        bulkp_cData = plot_data.show_data(indir, plot_type="bulkp_c",
            print_flag=False, save_flag=False, data_only=True)[1]
        fname = "bulkPotCathodeData.txt"
        np.savetxt(os.path.join(indir, fname), bulkp_cData,
                delimiter=",", header=header)

    return

if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise Exception("Need input directory name")
    indir = sys.argv[1]
    if not os.path.exists(os.path.join(os.getcwd(), indir)):
        raise Exception("Input directory doesn't exist")
    main(indir)
