import sys
import os

import numpy as np
import scipy.io as sio
import matplotlib as mpl
mpl.use("TkAgg")
import matplotlib.pyplot as plt
import matplotlib.animation as manim
import matplotlib.collections as mcollect

import mpet_params_IO
import elyte_CST

def show_data(indir, plot_type, print_flag, save_flag, data_only):
    pfx = 'mpet.'
    ttl_fmt = "% = {perc:2.1f}"
    # Read in the simulation results and calcuations data
    dataFileName = "output_data.mat"
    dataFile = os.path.join(indir, dataFileName)
    data = sio.loadmat(dataFile)
    try:
        data[pfx + 'current'][0][0]
    except KeyError:
        pfx = ''
    # Read in the parameters used to define the simulation
    paramFileName = "input_params_system.cfg"
    paramFile = os.path.join(indir, paramFileName)
    IO = mpet_params_IO.mpetIO()
    dD_s, ndD_s = IO.readDicts(os.path.join(indir, "input_dict_system"))
    # simulated (porous) electrodes
    Nvol = ndD_s["Nvol"]
    trodes = ndD_s["trodes"]
    dD_e = {}
    ndD_e = {}
    for trode in trodes:
        dD_e[trode], ndD_e[trode] = IO.readDicts(
                os.path.join(indir, "input_dict_{t}".format(t=trode)))
    # Pick out some useful constants/calculated values
    k = dD_s['k']                      # Boltzmann constant, J/(K Li)
    Tref = dD_s['Tref']                # Temp, K
    e = dD_s['e']                      # Charge of proton, C
    N_A = dD_s['N_A']                  # particles/mol
    F = dD_s['F']                      # C/mol
    td = dD_s["td"]
    Etheta = {"a" : 0.}
    for trode in trodes:
        Type = ndD_e[trode]['type']
        if Type not in ["diffn"]:
            Etheta[trode] = dD_e[trode]["Vstd"]
        else:
            Etheta[trode] = (k*Tref/e) * ndD_e[trode]["dphi_eq_ref"]
    Vstd = Etheta["c"] - Etheta["a"]
    Nvol = ndD_s["Nvol"]
    Npart = ndD_s["Npart"]
    psd_len = dD_s["psd_len"]
    # Discretization (and associated porosity)
    Lfac = 1e6
    Lunit = r"$\mu$m"
    dxc = ndD_s["L"]["c"]/Nvol["c"]
    dxvec = np.array(Nvol["c"] * [dxc])
    porosvec = np.array(Nvol["c"] * [ndD_s["poros"]["c"]])
    cellsvec = dxc*np.arange(Nvol["c"]) + dxc/2.
    if Nvol["s"]:
        dxs = ndD_s["L"]["s"]/Nvol["s"]
        dxvec_s = np.array(Nvol["s"] * [dxs])
        dxvec = np.hstack((dxvec_s, dxvec))
        poros_s = np.array(Nvol["s"] * [ndD_s["poros"]["s"]])
        porosvec = np.hstack((poros_s, porosvec))
        cellsvec += dD_s["L"]["s"] / dD_s["L"]["c"]
        cellsvec_s = dxs*np.arange(Nvol["s"]) + dxs/2.
        cellsvec = np.hstack((cellsvec_s, cellsvec))
    if "a" in trodes:
        dxa = ndD_s["L"]["a"]/Nvol["a"]
        dxvec_a = np.array(Nvol["a"] * [dxa])
        dxvec = np.hstack((dxvec_a, dxvec))
        poros_a = np.array(Nvol["a"] * [ndD_s["poros"]["a"]])
        porosvec = np.hstack((poros_a, porosvec))
        cellsvec += dD_s["L"]["a"] / dD_s["L"]["c"]
        cellsvec_a = dxa*np.arange(Nvol["a"]) + dxa/2.
        cellsvec = np.hstack((cellsvec_a, cellsvec))
    cellsvec *= dD_s["Lref"] * Lfac
    facesvec = np.insert(np.cumsum(dxvec), 0, 0.) * dD_s["Lref"] * Lfac
    # Extract the reported simulation times
    times = data[pfx + 'phi_applied_times'][0]
    numtimes = len(times)
    tmin = np.min(times)
    tmax = np.max(times)
    # Simulation type
    profileType = ndD_s['profileType']
    # Colors for plotting concentrations
    to_yellow = 0.3
    to_red = 0.7
    figsize = (6, 4.5)

    # Plot defaults
    axtickfsize = 18
    labelfsize = 20
    legfsize = labelfsize - 2
    txtfsize = labelfsize - 2
    lwidth = 3.
    markersize = 10
    mpl.rcParams['xtick.labelsize'] = axtickfsize
    mpl.rcParams['ytick.labelsize'] = axtickfsize
    mpl.rcParams['axes.labelsize'] = labelfsize
    mpl.rcParams['axes.labelsize'] = labelfsize
    mpl.rcParams['font.size'] = txtfsize
    mpl.rcParams['legend.fontsize'] = legfsize
    mpl.rcParams['lines.linewidth'] = lwidth
    mpl.rcParams['lines.markersize'] = markersize
    mpl.rcParams['lines.markeredgewidth'] = 0.1
#    mpl.rcParams['text.usetex'] = True

    # Print relevant simulation info
    if print_flag:
#        for i in trodes:
#            print "PSD_{l}:".format(l=l)
#            print psd_len[l].transpose()
#            print "Actual psd_mean [nm]:", np.mean(psd_len[l])
#            print "Actual psd_stddev [nm]:", np.std(psd_len[l])
        print "profileType:", profileType
        print "Cell structure:"
        print (("porous anode | " if Nvol["a"] else "flat anode | ")
                + ("sep | " if Nvol["s"] else "") + "porous cathode")
        if Nvol["a"]:
            print "capacity ratio cathode:anode, 'z':", ndD_s["z"]
        for trode in trodes:
            print "solidType_{t}:".format(t=trode), ndD_e[trode]['type']
            print "solidShape_{t}".format(t=trode), ndD_e[trode]['shape']
            print "rxnType_{t}:".format(t=trode), ndD_e[trode]['rxnType']
        if profileType == "CC":
            print "C_rate:", dD_s['Crate']
            print "current:", dD_s['currset'], "A/m^2"
        else: # CV
            print "Vset:", dD_s['Vset']
        print ("Specified psd_mean, c [{unit}]:".format(unit=Lunit),
                np.array(dD_s['psd_mean']["c"])*Lfac)
        print ("Specified psd_stddev, c [{unit}]:".format(unit=Lunit),
                np.array(dD_s['psd_stddev']["c"])*Lfac)
#        print "reg sln params:"
#        print ndD["Omga"]
        print "ndim B_c:", ndD_e["c"]["B"]
        if Nvol["s"]: print "Nvol_s:", Nvol["s"]
        print "Nvol_c:", Nvol["c"]
        if Nvol["a"]: print "Nvol_a:", Nvol["a"]
        print "Npart_c:", Npart["c"]
        if Nvol["a"]: print "Npart_a:", Npart["a"]
        print "Dp [m^2/s]:", dD_s['Dp']
        print "Dm [m^2/s]:", dD_s['Dm']
        print "Damb [m^2/s]:", dD_s['Damb']
        print "td [s]:", dD_s["td"]
        for trode in trodes:
            print "k0_{t} [A/m^2]:".format(t=trode), dD_e[trode]['k0']
            rxnType = ndD_e[trode]['rxnType']
            if rxnType == "BV":
                print "alpha_" + trode + ":", ndD_e[trode]['alpha']
            elif rxnType in ["Marcus", "MHC"]:
                print "lambda_" + trode + "/(kTref):", ndD_e[trode]["lambda"]
            if ndD_s['simBulkCond'][trode]:
                print (trode + " bulk conductivity loss: Yes -- " +
                        "dim_mcond [S/m]: " + str(dD_s['mcond'][trode]))
            else:
                print trode + " bulk conductivity loss: No"
            try:
                simSurfCond = ndD_e[trode]['simSurfCond']
                if simSurfCond:
                    print (trode + " surface conductivity loss: Yes -- " +
                            "dim_scond [S]: " + str(dD_e[trode]['scond']))
                else:
                    print trode + " surface conductivity loss: No"
            except:
                pass
#            if ndD['simSurfCond'][l]:
#                print (l + " surface conductivity loss: Yes -- " +
#                        "dim_scond [S]: " + str(dD['scond'][l]))
#            else:
#                print l + " surface conductivity loss: No"

    if plot_type in ["params"]:
        return ndD_s, dD_s, ndD_e, dD_e
    if plot_type in ["discData"]:
        return cellsvec/Lfac, facesvec/Lfac

    # Plot voltage profile
    if plot_type in ["v", "vt"]:
        voltage = (Vstd -
                (k*Tref/e)*data[pfx + 'phi_applied'][0])
        ffvec = data[pfx + 'ffrac_c'][0]
        fig, ax = plt.subplots()
        if plot_type == "v":
            if data_only:
                return ffvec, voltage
            ax.plot(ffvec, voltage)
            xmin = 0.
            xmax = 1.
            ax.set_xlim((xmin, xmax))
            ax.set_xlabel("Cathode Filling Fraction [dimensionless]")
        elif plot_type == "vt":
            if data_only:
                return times*td, voltage
            ax.plot(times*td, voltage)
            ax.set_xlabel("Time [s]")
        ax.set_ylabel("Voltage [V]")
        if save_flag:
            fig.savefig("mpet_v.png", bbox_inches="tight")
        return fig, ax

    # Plot surface conc.
    if plot_type in ["surf_c", "surf_a"]:
        l = plot_type[-1]
        if data_only:
            raise NotImplemented("no data-only output for surf")
        fig, ax = plt.subplots(Npart[l], Nvol[l], squeeze=False,
                sharey=True)
        str_base = pfx + "partTrode{l}vol{{j}}part{{i}}.".format(l=l) + "c"
        ylim = (0, 1.01)
        datax = times
        for i in range(Npart[l]):
            for j in range(Nvol[l]):
                sol_str = str_base.format(i=i, j=j)
                # Remove axis ticks
                ax[i, j].xaxis.set_major_locator(plt.NullLocator())
                datay = data[sol_str][:,-1]
                line, = ax[i, j].plot(times, datay)
        return fig, ax

    # Plot SoC profile
    if plot_type in ["soc_c", "soc_a"]:
        l = plot_type[-1]
        ffvec = data[pfx + 'ffrac_{l}'.format(l=l)][0]
        if data_only:
            return times*td, ffvec
        fig, ax = plt.subplots()
        print ffvec[-1]
        ax.plot(times*td, ffvec)
        xmin = np.min(ffvec)
        xmax = np.max(ffvec)
        ax.set_ylim((0, 1.05))
        ax.set_xlabel("Time [s]")
        ax.set_ylabel("Filling Fraciton [dimless]")
        if save_flag:
            fig.savefig("mpet_soc.png", bbox_inches="tight")
        return fig, ax

    # Check to make sure mass is conserved in elyte
    if plot_type == "elytecons":
        if data_only:
            raise NotImplemented("no data-only output for elytecons")
        fig, ax = plt.subplots()
        eps = 1e-2
        ymin = 1-eps
        ymax = 1+eps
#        ax.set_ylim((ymin, ymax))
        ax.set_ylabel('Avg. Concentration of electrolyte [nondim]')
        sep = pfx + 'c_lyte_s'
        anode = pfx + 'c_lyte_a'
        cath = pfx + 'c_lyte_c'
        ax.set_xlabel('Time [s]')
        cvec = data[cath]
        if Nvol["s"]:
            cvec_s = data[sep]
            cvec = np.hstack((cvec_s, cvec))
        if "a" in trodes:
            cvec_a = data[anode]
            cvec = np.hstack((cvec_a, cvec))
        cavg = np.sum(porosvec*dxvec*cvec,
                axis=1)/np.sum(porosvec*dxvec)
        np.set_printoptions(precision=8)
        print cavg
        ax.plot(times, cavg)
        return fig, ax

    # Plot current profile
    if plot_type == "curr":
        current = data[pfx + "current"][0] * 3600/td
        ffvec = data[pfx + 'ffrac_c'][0]
        if data_only:
            return times*td, current
        fig, ax = plt.subplots()
        ax.plot(times*td, current)
        xmin = np.min(ffvec)
        xmax = np.max(ffvec)
        ax.set_xlabel("Time [s]")
        ax.set_ylabel("Current [C-rate]")
        if save_flag:
            fig.savefig("mpet_current.png", bbox_inches="tight")
        return fig, ax

    # Plot electrolyte concentration or potential
    elif plot_type in ["elytec", "elytep", "elytecf", "elytepf",
            "elytei", "elyteif", "elytedivi", "elytedivif"]:
        fplot = (True if plot_type[-1] == "f" else False)
        t0ind = (0 if not fplot else -1)
        mpl.animation.Animation._blit_draw = _blit_draw
        datax = cellsvec
        c_sep, p_sep = pfx + 'c_lyte_s', pfx + 'phi_lyte_s'
        c_anode, p_anode = pfx + 'c_lyte_a', pfx + 'phi_lyte_a'
        c_cath, p_cath = pfx + 'c_lyte_c', pfx + 'phi_lyte_c'
        datay_c = data[c_cath]
        datay_p = data[p_cath]
        L_c = dD_s['L']["c"] * Lfac
        Ltot = L_c
        if Nvol["s"]:
            datay_s_c = data[c_sep]
            datay_s_p = data[p_sep]
            datay_c = np.hstack((datay_s_c, datay_c))
            datay_p = np.hstack((datay_s_p, datay_p))
            L_s = dD_s['L']["s"] * Lfac
            Ltot += L_s
        else:
            L_s = 0
        if "a" in trodes:
            datay_a_c = data[c_anode]
            datay_a_p = data[p_anode]
            datay_c = np.hstack((datay_a_c, datay_c))
            datay_p = np.hstack((datay_a_p, datay_p))
            L_a = dD_s['L']["a"] * Lfac
            Ltot += L_a
        else:
            L_a = 0
        xmin = 0
        xmax = Ltot
        if plot_type in ["elytec", "elytecf"]:
            ymin = 0
            ymax = 2.2
            ylbl = 'Concentration of electrolyte [M]'
            datay = datay_c * dD_s["cref"] / 1000.
        elif plot_type in ["elytep", "elytepf"]:
            ymin = -50
            ymax = 50
            ylbl = 'Potential of electrolyte [V]'
            datay = datay_p * k*Tref/e
        elif plot_type in ["elytei", "elyteif", "elytedivi",
            "elytedivif"]:
            dxd1 = (dxvec[0:-1] + dxvec[1:]) / 2.
            dxd2 = dxvec
            c_edges = (2*datay_c[:, 1:]*datay_c[:, :-1])/(datay_c[:, 1:] +
                    datay_c[:, :-1] + 1e-20)
            porosTmp = porosvec**(1.5)
            poros_edges = (2*porosTmp[1:]*porosTmp[:-1])/(porosTmp[1:]
                    + porosTmp[:-1] + 1e-20)
            if ndD_s["elyteModelType"] == "SM":
                D_lyte, kappa_lyte, thermFac_lyte, tp0, Dref = (
                        elyte_CST.getProps(ndD_s["SMset"]))
                i_edges = -poros_edges*kappa_lyte(c_edges) * (
                        np.diff(datay_p, axis=1)/dxd1 -
                        (ndD_s["nup"]+ndD_s["num"])/ndD_s["nup"] *
                        (1 - tp0(c_edges)) *
                        thermFac_lyte(c_edges) *
                        np.diff(np.log(datay_c), axis=1)/dxd1
                        )
            elif ndD_s["elyteModelType"] == "dilute":
                Dp, Dm = ndD_s["Dp"], ndD_s["Dm"]
                zp, zm = ndD_s["zp"], ndD_s["zm"]
                i_edges = (-(Dp - Dm)*np.diff(datay_c, axis=1)/dxd1
                        - (zp*Dm - zm*Dm)*c_edges*np.diff(datay_p, axis=1)/dxd1
                        )
            i_CCs = np.zeros((numtimes, 1))
            i_edges = np.hstack((i_CCs, i_edges, i_CCs))
            if plot_type in ["elytei", "elyteif"]:
                ylbl = r'Current density of electrolyte [A/m$^2$]'
                datax = facesvec
                datay = i_edges * (F*dD_s["cref"]*dD_s["Dref"]/dD_s["Lref"])
            elif plot_type in ["elytedivi", "elytedivif"]:
                ylbl = r'Divergence of Current density of electrolyte [A/m$^3$]'
                datax = cellsvec
                datay = np.diff(i_edges, axis=1) / dxd2
                datay *= (F*dD_s["cref"]*dD_s["Dref"]/dD_s["Lref"]**2)
        if fplot:
            datay = datay[t0ind]
        if data_only:
            return datax, datay
        dataMin, dataMax = np.min(datay), np.max(datay)
        dataRange = dataMax - dataMin
        ymin = dataMin - 0.05*dataRange
        ymax = dataMax + 0.05*dataRange
        fig, ax = plt.subplots()
        ax.set_xlabel('Battery Position [{unit}]'.format(unit=Lunit))
        ax.set_ylabel(ylbl)
        ttl = ax.text(0.5, 1.05, ttl_fmt.format(perc=0),
                transform = ax.transAxes, verticalalignment="center",
                horizontalalignment="center")
        ax.set_ylim((ymin, ymax))
        ax.set_xlim((xmin, xmax))
        # returns tuble of line objects, thus comma
        line1, = ax.plot(datax, datay[t0ind, :], '-')
        ax.axvline(x=L_a, linestyle='--', color='g')
        ax.axvline(x=(L_a+L_s), linestyle='--', color='g')
        if fplot:
            return fig, ax
        def init():
            line1.set_ydata(np.ma.array(datax, mask=True))
            ttl.set_text('')
            return line1, ttl
        def animate(tind):
            line1.set_ydata(datay[tind])
            t_current = times[tind]
            tfrac = (t_current - tmin)/(tmax - tmin) * 100
            ttl.set_text(ttl_fmt.format(perc=tfrac))
            return line1, ttl

    # Plot all solid concentrations or potentials
    elif plot_type in ["csld_c", "csld_a", "phisld_a", "phisld_c"]:
        t0ind = 0
        l = plot_type[-1]
        if data_only:
            raise NotImplemented("no data-only output for csld/phisld")
        fig, ax = plt.subplots(Npart[l], Nvol[l], squeeze=False,
                sharey=True)
        sol = np.empty((Npart[l], Nvol[l]), dtype=object)
        sol1 = np.empty((Npart[l], Nvol[l]), dtype=object)
        sol2 = np.empty((Npart[l], Nvol[l]), dtype=object)
        lens = np.zeros((Npart[l], Nvol[l]))
        lines = np.empty((Npart[l], Nvol[l]), dtype=object)
        lines1 = np.empty((Npart[l], Nvol[l]), dtype=object)
        lines2 = np.empty((Npart[l], Nvol[l]), dtype=object)
        partStr = "partTrode{l}vol{j}part{i}."
        type2c = False
        if plot_type in ["csld_a", "csld_col_a", "csld_c", "csld_col_c"]:
            if ndD_e[l]["type"] in ndD_s["1varTypes"]:
                str_base = pfx + partStr + "c"
            elif ndD_e[l]["type"] in ndD_s["2varTypes"]:
                type2c = True
                str1_base = pfx + partStr + "c1"
                str2_base = pfx + partStr + "c2"
            ylim = (0, 1.01)
        elif plot_type in ["phisld_a", "phisld_c"]:
            str_base = pfx + partStr + "phi"
            ylim = (-10, 20)
        for i in range(Npart[l]):
            for j in range(Nvol[l]):
                lens[i, j] = psd_len[l][j, i]
                if type2c:
                    sol1[i, j] = str1_base.format(l=l, i=i, j=j)
                    sol2[i, j] = str2_base.format(l=l, i=i, j=j)
                    datay1 = data[sol1[i, j]][t0ind]
                    datay2 = data[sol2[i, j]][t0ind]
                    numy = len(datay1)
                    datax = np.linspace(0, lens[i, j], numy)
                    line1, = ax[i, j].plot(datax, datay1)
                    line2, = ax[i, j].plot(datax, datay2)
                    lines1[i, j] = line1
                    lines2[i, j] = line2
                else:
                    sol[i, j] = str_base.format(l=l, i=i, j=j)
                    datay = data[sol[i, j]][t0ind]
                    numy = len(datay)
                    datax = np.linspace(0, lens[i, j], numy)
                    line, = ax[i, j].plot(datax, datay)
                    lines[i, j] = line
                # Remove axis ticks
                ax[i, j].set_ylim(ylim)
                ax[i, j].set_xlim((0, lens[i, j]))
        def init():
            for i in range(Npart[l]):
                for j in range(Nvol[l]):
                    if type2c:
                        numy = len(data[sol1[i, j]][t0ind])
                        maskTmp = np.zeros(numy)
                        lines1[i, j].set_ydata(np.ma.array(maskTmp, mask=True))
                        lines2[i, j].set_ydata(np.ma.array(maskTmp, mask=True))
                    else:
                        numy = len(data[sol[i, j]][t0ind])
                        maskTmp = np.zeros(numy)
                        lines[i, j].set_ydata(np.ma.array(maskTmp, mask=True))
            if type2c:
                return tuple(np.vstack((lines1, lines2)).reshape(-1))
            else:
                return tuple(lines.reshape(-1))
        def animate(tind):
            for i in range(Npart[l]):
                for j in range(Nvol[l]):
                    if type2c:
                        datay1 = data[sol1[i, j]][tind]
                        datay2 = data[sol2[i, j]][tind]
                        lines1[i, j].set_ydata(datay1)
                        lines2[i, j].set_ydata(datay2)
                    else:
                        datay = data[sol[i, j]][tind]
                        lines[i, j].set_ydata(datay)
            if type2c:
                return tuple(np.vstack((lines1, lines2)).reshape(-1))
            else:
                return tuple(lines.reshape(-1))

    # Plot average solid concentrations
    elif plot_type in ["cbar_c", "cbar_a", "cbar_full"]:
        if plot_type[-1] == "l":
            lvec = ["a", "c"]
        elif plot_type[-1] == "a":
            lvec = ["a"]
        else:
            lvec = ["c"]
        dataCbar = {}
        for trode in trodes:
            dataCbar[trode] = np.zeros((numtimes, Nvol[trode], Npart[trode]))
            for tInd in range(numtimes):
                for i in range(Nvol[trode]):
                    for j in range(Npart[trode]):
                        dataStr = (pfx +
                                'partTrode{t}vol{i}part{j}.'.format(t=trode,i=i,j=j)
                                + 'cbar')
                        dataCbar[trode][tInd, i, j] = (
                                data[dataStr][0][tInd])
        if data_only:
            return dataCbar
        # Set up colors.
        # Define if you want smooth or discrete color changes
        # Option: "smooth" or "discrete"
        color_changes = "discrete"
#        color_changes = "smooth"
        # Discrete color changes:
        if color_changes == "discrete":
            # Make a discrete colormap that goes from green to yellow
            # to red instantaneously
            cdict = {
                    "red" : [(0.0, 0.0, 0.0),
                             (to_yellow, 0.0, 1.0),
                             (1.0, 1.0, 1.0)],
                    "green" : [(0.0, 0.502, 0.502),
                               (to_yellow, 0.502, 1.0),
                               (to_red, 1.0, 0.0),
                               (1.0, 0.0, 0.0)],
                    "blue" : [(0.0, 0.0, 0.0),
                              (1.0, 0.0, 0.0)]
                    }
            cmap = mpl.colors.LinearSegmentedColormap(
                    "discrete", cdict)
        # Smooth colormap changes:
        if color_changes == "smooth":
#            cmap = mpl.cm.RdYlGn_r # A default green-yellow-red map
            # generated with colormap.org
            cmaps = np.load("colormaps_custom.npz")
            cmap_data = cmaps["GnYlRd_3"]
            cmap = mpl.colors.ListedColormap(cmap_data/255.)

        # Implement hack to be able to animate title
        mpl.animation.Animation._blit_draw = _blit_draw
        size_frac_min = 0.10
        fig, axs = plt.subplots(1, len(lvec), squeeze=False)
        ttlx = 0.5 if len(lvec) < 2 else 1.1
        ttl = axs[0,0].text(ttlx, 1.05, ttl_fmt.format(perc=0),
                transform = axs[0,0].transAxes, verticalalignment="center",
                horizontalalignment="center")
        collection = np.empty(len(lvec), dtype=object)
        for indx, l in enumerate(lvec):
            ax = axs[0, indx]
            # Get particle sizes (and max size) (length-based)
            lens = psd_len[l]
            len_max = np.max(lens)
            len_min = np.min(lens)
            ax.patch.set_facecolor('white')
            # Don't stretch axes to fit figure -- keep 1:1 x:y ratio.
            ax.set_aspect('equal', 'box')
            # Don't show axis ticks
            ax.xaxis.set_major_locator(plt.NullLocator())
            ax.yaxis.set_major_locator(plt.NullLocator())
            ax.set_xlim(0, 1.)
            ax.set_ylim(0, float(Npart[l])/Nvol[l])
            # Label parts of the figure
#            ylft = ax.text(-0.07, 0.5, "Separator",
#                    transform=ax.transAxes, rotation=90,
#                    verticalalignment="center",
#                    horizontalalignment="center")
#            yrht = ax.text(1.09, 0.5, "Current Collector",
#                    transform=ax.transAxes, rotation=90,
#                    verticalalignment="center",
#                    horizontalalignment="center")
#            xbtm = ax.text(.50, -0.05, "Electrode Depth -->",
#                    transform=ax.transAxes, rotation=0,
#                    verticalalignment="center",
#                    horizontalalignment="center")
            # Geometric parameters for placing the rectangles on the axes
            spacing = 1.0 / Nvol[l]
            size_fracs = 0.4*np.ones((Nvol[l], Npart[l]))
            if len_max != len_min:
                size_fracs = (lens - len_min)/(len_max - len_min)
            sizes = (size_fracs * (1 - size_frac_min) + size_frac_min) / Nvol[l]
            # Create rectangle "patches" to add to figure axes.
            rects = np.empty((Nvol[l], Npart[l]), dtype=object)
            color = 'green' # value is irrelevant -- it will be animated
            for (i, j), c in np.ndenumerate(sizes):
                size = sizes[i, j]
                center = np.array([spacing*(i + 0.5), spacing*(j + 0.5)])
                bottom_left = center - size / 2
                rects[i, j] = plt.Rectangle(bottom_left,
                        size, size, color=color)
            # Create a group of rectange "patches" from the rects array
            collection[indx] = mcollect.PatchCollection(rects.reshape(-1))
            # Put them on the axes
            ax.add_collection(collection[indx])
        # Have a "background" image of rectanges representing the
        # initial state of the system.
        def init():
            for indx, l in enumerate(lvec):
                cbar_mat = dataCbar[l][0, :, :]
                colors = cmap(cbar_mat.reshape(-1))
                collection[indx].set_color(colors)
                ttl.set_text('')
            out = [collection[i] for i in range(len(collection))]
            out.append(ttl)
            out = tuple(out)
            return out
        def animate(tind):
            for indx, l in enumerate(lvec):
                cbar_mat = dataCbar[l][tind, :, :]
                colors = cmap(cbar_mat.reshape(-1))
                collection[indx].set_color(colors)
            t_current = times[tind]
            tfrac = (t_current - tmin)/(tmax - tmin) * 100
            ttl.set_text(ttl_fmt.format(perc=tfrac))
            out = [collection[i] for i in range(len(collection))]
            out.append(ttl)
            out = tuple(out)
            return out

    # Plot cathode potential
    elif plot_type in ["bulkp_c", "bulkp_a"]:
        l = plot_type[-1]
        t0ind = 0
        mpl.animation.Animation._blit_draw = _blit_draw
        fig, ax = plt.subplots()
        ax.set_xlabel('Position in electrode [{unit}]'.format(unit=Lunit))
        ax.set_ylabel('Potential of cathode [nondim]')
        ttl = ax.text(0.5, 1.05, ttl_fmt.format(perc=0),
                transform = ax.transAxes, verticalalignment="center",
                horizontalalignment="center")
        bulkp = pfx + 'phi_bulk_{l}'.format(l=l)
        Ltrode = dD_s['L'][l] * Lfac
        datay = data[bulkp]
        ymin = np.min(datay) - 0.2
        ymax = np.max(datay) + 0.2
        if l == "a":
            datax = cellsvec[:Nvol["a"]]
        elif l == "c":
            datax = cellsvec[-Nvol["c"]:]
        if data_only:
            return datax, datay
        # returns tuble of line objects, thus comma
        line1, = ax.plot(datax, datay[t0ind])
        def init():
            line1.set_ydata(np.ma.array(datax, mask=True))
            ttl.set_text('')
            return line1, ttl
        def animate(tind):
            line1.set_ydata(datay[tind])
            t_current = times[tind]
            tfrac = (t_current - tmin)/(tmax - tmin) * 100
            ttl.set_text(ttl_fmt.format(perc=tfrac))
            return line1, ttl

    else:
        raise Exception("Unexpected plot type argument. " +
                "Try 'v', 'curr', 'elytec', 'elytep', " +
                "'soc_a/c', " +
                "'cbar_a/c', 'csld_a/c', 'phisld_a/c', " +
                "'bulkp_a/c', 'surf_a/c'.")

    ani = manim.FuncAnimation(fig, animate, frames=numtimes,
            interval=50, blit=True, repeat=False, init_func=init)
    if save_flag:
        ani.save("mpet_{type}.mp4".format(type=plot_type),
                fps=25, bitrate=5500,
#                writer='alzkes',
#                savefig_kwargs={'bbox_inches' : 'tight'},
                )
#                extra_args=['-vcodec', 'libx264'])

    return fig, ax, ani

# This is a block of code which messes with some mpl internals
# to allow for animation of a title. See
# http://stackoverflow.com/questions/17558096/animated-title-in-mpl
def _blit_draw(self, artists, bg_cache):
    # Handles blitted drawing, which renders only the artists given instead
    # of the entire figure.
    updated_ax = []
    for a in artists:
        # If we haven't cached the background for this axes object, do
        # so now. This might not always be reliable, but it's an attempt
        # to automate the process.
        if a.axes not in bg_cache:
            # bg_cache[a.axes] = a.figure.canvas.copy_from_bbox(a.axes.bbox)
            # change here
            bg_cache[a.axes] = a.figure.canvas.copy_from_bbox(a.axes.figure.bbox)
        a.axes.draw_artist(a)
        updated_ax.append(a.axes)

    # After rendering all the needed artists, blit each axes individually.
    for ax in set(updated_ax):
        # and here
        # ax.figure.canvas.blit(ax.bbox)
        ax.figure.canvas.blit(ax.figure.bbox)

if __name__ == "__main__":
    # Get input file from script parameters
    if len(sys.argv) < 2:
        raise Exception("Need input data directory name")
    indir = sys.argv[1]
    if not os.path.exists(os.path.join(os.getcwd(), indir)):
        raise Exception("Input file doesn't exist")
    # Get plot type from script parameters
    plots = []
    if len(sys.argv) > 2:
        plots.append(sys.argv[2])
    else:
        plots.append("v")
    # Save the plot instead of showing on screen?
    # Get from script parameters
    save_flag = False
    print_flag = True
    data_only = False
    if len(sys.argv) > 3:
        if sys.argv[3] == "save":
            save_flag = True
        else:
            for i in range(3, len(sys.argv)):
                plots.append(sys.argv[i])
    out = []
    for plot_type in plots:
        out.append(show_data(indir, plot_type, print_flag, save_flag,
            data_only))
    plt.show()
