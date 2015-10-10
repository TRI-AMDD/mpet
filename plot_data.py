import sys
import os
import time

import numpy as np
import scipy.io as sio
import matplotlib as mpl
mpl.use("TkAgg")
import matplotlib.pyplot as plt
import matplotlib.animation as manim
import matplotlib.collections as mcollect

import mpet_params_IO

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
    paramFileName = "input_params.cfg"
    paramFile = os.path.join(indir, paramFileName)
    IO = mpet_params_IO.mpetIO()
    dD, ndD = IO.readDicts(os.path.join(indir, "input_dict"))
    # simulated (porous) electrodes
    Nvol = ndD["Nvol"]
    trodes = ndD["trodes"]
    # Pick out some useful constants/calculated values
    k = dD['k']                      # Boltzmann constant, J/(K Li)
    Tref = dD['Tref']                # Temp, K
    e = dD['e']                      # Charge of proton, C
    td = dD["td"]
    Vstd = dD["Vstd"]
    for l in trodes:
        # Replace the standard potential if a fit voltage curve was used.
        # Use the value that the simulation used in initialization.
        if ndD['delPhiEqFit'][l]:
            Vstd[l] = ndD["dphi_eq_ref"][l]*(k*Tref/e)
    Nvol = ndD["Nvol"]
    Npart = ndD["Npart"]
    solidType = ndD["solidType"]
    solidShape = ndD["solidShape"]
    rxnType = ndD["rxnType"]
    psd_len = dD["psd_len"]
    # Discretization (and associated porosity)
    Lfac = 1e6
    Lunit = r"$\mu$m"
    dxc = ndD["L"]["c"]/Nvol["c"]
    dxvec = np.array(Nvol["c"] * [dxc])
    porosvec = np.array(Nvol["c"] * [ndD["poros"]["c"]])
#    L_c = dD['L']["c"] * Lfac
    cellsvec = dxc*np.arange(Nvol["c"]) + dxc/2.
    if Nvol["s"]:
        dxs = ndD["L"]["s"]/Nvol["s"]
        dxvec_s = np.array(Nvol["s"] * [dxs])
        dxvec = np.hstack((dxvec_s, dxvec))
        poros_s = np.array(Nvol["s"] * [ndD["poros"]["s"]])
        porosvec = np.hstack((poros_s, porosvec))
#        L_s = dD['L']["s"] * Lfac
        cellsvec += dD["L"]["s"] / dD["L"]["c"]
        cellsvec_s = dxs*np.arange(Nvol["s"]) + dxs/2.
        cellsvec = np.hstack((cellsvec_s, cellsvec))
    if "a" in trodes:
        dxa = ndD["L"]["a"]/Nvol["a"]
        dxvec_a = np.array(Nvol["a"] * [dxa])
        dxvec = np.hstack((dxvec_a, dxvec))
        poros_a = np.array(Nvol["a"] * [ndD["poros"]["a"]])
        porosvec = np.hstack((poros_a, porosvec))
#        L_a = dD['L']["a"] * Lfac
        cellsvec += dD["L"]["a"] / dD["L"]["c"]
        cellsvec_a = dxa*np.arange(Nvol["a"]) + dxa/2.
        cellsvec = np.hstack((cellsvec_a, cellsvec))
    cellsvec *= dD["L"]["c"] * Lfac
    # Extract the reported simulation times
    times = data[pfx + 'phi_applied_times'][0]
    numtimes = len(times)
    tmin = np.min(times)
    tmax = np.max(times)
    # Simulation type
    profileType = ndD['profileType']
    # Colors for plotting concentrations
    to_yellow = 0.3
    to_red = 0.7

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
            print "capacity ratio cathode:anode, 'z':", ndD["z"]
        print "solidType:", solidType
        print "solidShape", solidShape
        print "rxnType:", rxnType
        if profileType == "CC":
            print "C_rate:", dD['Crate']
        else: # CV
            print "Vset:", dD['Vset']
        print ("Specified psd_mean, c [{unit}]:".format(unit=Lunit),
                np.array(dD['psd_mean']["c"])*Lfac)
        print ("Specified psd_stddev, c [{unit}]:".format(unit=Lunit),
                np.array(dD['psd_stddev']["c"])*Lfac)
#        print "reg sln params:"
#        print ndD["Omga"]
        if Nvol["s"]: print "Nvol_s:", Nvol["s"]
        print "Nvol_c:", Nvol["c"]
        if Nvol["a"]: print "Nvol_a:", Nvol["a"]
        print "Npart_c:", Npart["c"]
        if Nvol["a"]: print "Npart_a:", Npart["a"]
        print "Dp [m^2/s]:", dD['Dp']
        print "Dm [m^2/s]:", dD['Dm']
        print "Damb [m^2/s]:", dD['Damb']
        print "td [s]:", dD["td"]
        print "k0 [A/m^2]:", dD['k0']
        for l in trodes:
            if rxnType[l] == "BV":
                print "alpha_" + l + ":", ndD['alpha'][l]
            elif rxnType[l] in ["Marcus", "MHC"]:
                print "lambda_" + l + "/(kTref):", ndD["lambda"][l]
            if ndD['simBulkCond'][l]:
                print (l + " bulk conductivity loss: Yes -- " +
                        "dim_mcond [S/m]: " + str(dD['mcond'][l]))
            else:
                print l + " bulk conductivity loss: No"
            if ndD['simSurfCond'][l]:
                print (l + " surface conductivity loss: Yes -- " +
                        "dim_scond [S]: " + str(dD['scond'][l]))
            else:
                print l + " surface conductivity loss: No"

    # Plot voltage profile
    if plot_type in ["v", "vt"]:
        voltage = ((Vstd["c"] - Vstd["a"]) -
                (k*Tref/e)*data[pfx + 'phi_applied'][0])
        ffvec = data[pfx + 'ffrac_c'][0]
        if data_only:
            return ffvec, voltage
        fig, ax = plt.subplots()
        if plot_type == "v":
            ax.plot(ffvec, voltage)
            xmin = 0.
            xmax = 1.
            ax.set_xlim((xmin, xmax))
            ax.set_xlabel("Cathode Filling Fraction [dimensionless]")
        elif plot_type == "vt":
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
        str_base = pfx + "c_sld_trode{l}vol{{j}}part{{i}}".format(l=l)
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

    # Plot misc stuff about reactions
    if plot_type in ["rxnp_c", "rxnp_a"]:
        l = plot_type[-1]
        if data_only:
            raise NotImplemented("no data-only output for rxnp")
        fig, ax = plt.subplots(Npart[l], Nvol[l], squeeze=False,
                sharey=True)
        k0 = dD['k0'][l]
        sol_c_str_base = pfx + "c_sld_trode{l}vol{{j}}part{{i}}".format(l=l)
        sol_p_str = pfx + "phi_{l}".format(l=l)
        lyte_c_str = pfx + "c_lyte_{l}".format(l=l)
        lyte_p_str = pfx + "phi_lyte_{l}".format(l=l)
        ylim = (0, 1.01)
        ffvec = data[pfx + 'ffrac_{l}'.format(l=l)][0]
        datax = times
#        datax = ffvec
        for i in range(Npart[l]):
            for j in range(Nvol[l]):
                sol_str = sol_c_str_base.format(i=i, j=j)
                # Remove axis ticks
#                ax[i, j].xaxis.set_major_locator(plt.NullLocator())
                csld = data[sol_str]
                c_lyte = data[lyte_c_str][:,j]
                phi_lyte = data[lyte_p_str][:,j]
                phi_m = data[sol_p_str][:,j]
#                csld = csld[:, 50]
                Omga = data[pfx + "Omga_{l}".format(l=l)][0][j][i]
                def mu_reg_sln(csld, Omga):
                    return Omga*(1-2*csld) + np.log(csld/(1-csld))
#                # homog particles
#                mu_R = Omga*(1-2*csld) + 1*np.log(csld/(1-csld))
#                # ACR particles
#                sld_pt_indx = 50
#                Nij = csld.shape[1]
#                cbar = np.tile(np.sum(csld, axis=1)/Nij, [Nij, 1]).transpose()
#                cstmp = np.zeros((len(datax), Nij + 2))
#                cstmp[:, 1:-1] = csld
#                cstmp[0] = D['cwet_ac'][l]
#                cstmp[-1] = D['cwet_ac'][l]
#                dxs = 1./Nij
#                curv = np.diff(cstmp,2)/(dxs**2)
#                kappa = data[pfx + "kappa_{l}".format(l=l)][i, j]
#                B = data[pfx + "B_ac"][0][l]
#                mu_R = ( mu_reg_sln(csld, Omga) - kappa*curv
#                    + B*(csld - cbar) )
#                mu_R = mu_R[:, sld_pt_indx]
#                act_R = np.exp(mu_R)
#                eta = (mu_R + phi_m) - (mu_O + phi_lyte)
#                BVgamma_ts = 1./(1-csld[:, sld_pt_indx])
                # diffn with delPhiEqFit
                csurf = csld[:,-1]
                import delta_phi_fits
                dphi_eq_ref = data[pfx + "dphi_eq_ref_ac"][0][l]
                fits = delta_phi_fits.DPhiFits(D)
                phifunc = fits.materialData[D["material_ac"][l]]
                dphi_eq = phifunc(csurf, dphi_eq_ref)
                eta = phi_m - phi_lyte - dphi_eq
                act_R = csurf
                BVgamma_ts = 1./(1-csurf)
                # General
                act_O = c_lyte
                mu_O = np.log(act_O)
#                eta = np.linspace(-20, 20, 70)
                BValpha = ndD['alpha'][l]
                BVecd = ( k0 * act_O**(1-BValpha)
                    * act_R**(BValpha) / BVgamma_ts )
                BVrate = ( BVecd *
                    (np.exp(-BValpha*eta/1) - np.exp((1-BValpha)*eta/1)) )
#                Malpha = 0.5*(1 + (1/lmbda) * np.log(c_lyte/csld))
#                Mecd = ( k0 *
#                    c_lyte**((3-2*Malpha)/4.) *
#                    csld**((1+2*Malpha)/4.) )
#                Meta2 = np.exp(-eta**2/(4.*1*lmbda))
#                Mrate = ( Mecd * np.exp(-eta**2/(4.*1*lmbda)) *
#                    (np.exp(-Malpha*eta/1) - np.exp((1-Malpha)*eta/1)) )
                line, = ax[i, j].plot(datax, BVecd)
        if save_flag:
            fig.savefig("Rxn_out.png", bbox_inches="tight")
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
#        print times*td, current
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
    elif plot_type == "elytec" or plot_type == "elytep":
        mpl.animation.Animation._blit_draw = _blit_draw
        if data_only:
            raise NotImplemented("no data-only output for elytec/p")
        fig, ax = plt.subplots()
        if plot_type == "elytec":
            ymin = 0
            ymax = 2.2
            ax.set_ylabel('Concentration of electrolyte [nondim]')
            sep = pfx + 'c_lyte_s'
            anode = pfx + 'c_lyte_a'
            cath = pfx + 'c_lyte_c'
        elif plot_type == "elytep":
            ymin = -50
            ymax = 50
            ax.set_ylabel('Potential of electrolyte [nondim]')
            sep = pfx + 'phi_lyte_s'
            anode = pfx + 'phi_lyte_a'
            cath = pfx + 'phi_lyte_c'
        ax.set_xlabel('Battery Position [{unit}]'.format(unit=Lunit))
        ttl = ax.text(0.5, 1.05, ttl_fmt.format(perc=0),
                transform = ax.transAxes, verticalalignment="center",
                horizontalalignment="center")
        datay = data[cath][0]
        L_c = dD['L']["c"] * Lfac
        Ltot = L_c
        if Nvol["s"]:
            datay_s = data[sep][0]
            datay = np.hstack((datay_s, datay))
            L_s = dD['L']["s"] * Lfac
            Ltot += L_s
        else:
            L_s = 0
        if "a" in trodes:
            datay_a = data[anode][0]
            datay = np.hstack((datay_a, datay))
            L_a = dD['L']["a"] * Lfac
            Ltot += L_a
        else:
            L_a = 0
        numy = len(datay)
        xmin = 0
        xmax = Ltot
        datax = cellsvec
        ax.set_ylim((ymin, ymax))
        ax.set_xlim((xmin, xmax))
        # returns tuble of line objects, thus comma
        line1, = ax.plot(datax, datay, '-')
        ax.axvline(x=L_a, linestyle='--', color='g')
        ax.axvline(x=(L_a+L_s), linestyle='--', color='g')
        def init():
            line1.set_ydata(np.ma.array(datax, mask=True))
            ttl.set_text('')
            return line1, ttl
        def animate(tind):
            datay = data[cath][tind]
            if Nvol["s"]:
                datay_s = data[sep][tind]
                datay = np.hstack((datay_s, datay))
            if "a" in trodes:
                datay_a = data[anode][tind]
                datay = np.hstack((datay_a, datay))
            line1.set_ydata(datay)
            t_current = times[tind]
            tfrac = (t_current - tmin)/(tmax - tmin) * 100
            ttl.set_text(ttl_fmt.format(perc=tfrac))
            return line1, ttl

    # Plot all solid concentrations or potentials
    elif plot_type in ["csld_c", "csld_a", "phisld_a", "phisld_c",
            "csld_col_c", "csld_col_a"]:
        l = plot_type[-1]
        if data_only:
            raise NotImplemented("no data-only output for csld/phisld")
        fig, ax = plt.subplots(Npart[l], Nvol[l], squeeze=False,
                sharey=True)
        sol = np.empty((Npart[l], Nvol[l]), dtype=object)
        lens = np.zeros((Npart[l], Nvol[l]))
        lines = np.empty((Npart[l], Nvol[l]), dtype=object)
        fills = np.empty((Npart[l], Nvol[l], 3), dtype=object)
        if plot_type in ["csld_a", "csld_col_a", "csld_c", "csld_col_c"]:
            str_base = pfx + "c_sld_trode{l}vol{j}part{i}"
            ylim = (0, 1.01)
#        elif plot_type in ["csld_c", "csld_col_c"]:
#            str_base = pfx + "c_sld_trode1vol{j}part{i}"
#            ylim = (0, 1.01)
        elif plot_type in ["phisld_a", "phisld_c"]: # plot_type == "phisld"
            str_base = pfx + "p_sld_trode{l}vol{j}part{i}"
            ylim = (-10, 20)
#        elif plot_type == "phisld_c":
#            str_base = pfx + "p_sld_trode1vol{j}part{i}"
#            ylim = (-10, 20)
        for i in range(Npart[l]):
            for j in range(Nvol[l]):
                sol[i, j] = str_base.format(l=l, i=i, j=j)
                lens[i, j] = psd_len[l][j, i]
                # Remove axis ticks
                ax[i, j].xaxis.set_major_locator(plt.NullLocator())
#                ax[i, j].yaxis.set_major_locator(plt.NullLocator())
                datay = data[sol[i, j]][0]
                numy = len(datay)
                datax = np.linspace(0, lens[i, j], numy)
                ax[i, j].set_ylim(ylim)
                ax[i, j].set_xlim((0, lens[i, j]))
                line, = ax[i, j].plot(datax, datay)
                lines[i, j] = line
                if plot_type in ["csld_col_c", "csld_col_a"]:
                    fill1 = ax[i, j].fill_between(datax, ylim[0],
                            ylim[1], facecolor='red', alpha=0.9,
                            where=datay>to_red)
                    fill2 = ax[i, j].fill_between(datax, ylim[0],
                            ylim[1], facecolor='yellow', alpha=0.9,
                            where=((datay<to_red) & (datay>to_yellow)))
                    fill3 = ax[i, j].fill_between(datax, ylim[0],
                            ylim[1], facecolor='green', alpha=0.9,
                            where=datay<to_yellow)
                    fills[i, j, :] = [fill1, fill2, fill3]
        def init():
            for i in range(Npart[l]):
                for j in range(Nvol[l]):
                    datax = np.zeros(data[sol[i, j]][0].shape)
                    lines[i, j].set_ydata(np.ma.array(datax, mask=True))
                    if plot_type in ["csld_col_c", "csld_col_a"]:
                        fill1 = ax[i, j].fill_between(datax, ylim[0],
                                ylim[1], facecolor='red', alpha=0.0)
                        fill2 = ax[i, j].fill_between(datax, ylim[0],
                                ylim[1], facecolor='yellow', alpha=0.0)
                        fill3 = ax[i, j].fill_between(datax, ylim[0],
                                ylim[1], facecolor='green', alpha=0.0)
                        fills[i, j, :] = [fill1, fill2, fill3]
#            if plot_type == "csld_col":
            if plot_type in ["csld_col_c", "csld_col_a"]:
#                collection = mcollect.PatchCollection(fills.reshape(-1))
#                return tuple(collection)
                return tuple(fills.reshape(-1))
            else:
                return tuple(lines.reshape(-1))
        def animate(tind):
            for i in range(Npart[l]):
                for j in range(Nvol[l]):
                    datay = data[sol[i, j]][tind]
                    lines[i, j].set_ydata(datay)
                    datax = lines[i, j].get_xdata()
                    if plot_type in ["csld_col_c", "csld_col_a"]:
                        fill1 = ax[i, j].fill_between(datax, ylim[0],
                                ylim[1], facecolor='red', alpha=0.9,
                                where=datay>to_red)
                        fill2 = ax[i, j].fill_between(datax, ylim[0],
                                ylim[1], facecolor='yellow', alpha=0.9,
                                where=((datay<to_red) & (datay>to_yellow)))
                        fill3 = ax[i, j].fill_between(datax, ylim[0],
                                ylim[1], facecolor='green', alpha=0.9,
                                where=datay<to_yellow)
                        fills[i, j, :] = [fill1, fill2, fill3]
#            if plot_type == "csld_col":
            if plot_type in ["csld_col_c", "csld_col_a"]:
#                collection = mcollect.PatchCollection(fills.reshape(-1))
#                return tuple(collection)
                return tuple(fills.reshape(-1))
            else:
                return tuple(lines.reshape(-1))

    # Plot average solid concentrations
    elif plot_type in ["cbar_c", "cbar_a", "cbar_full"]:
        if data_only:
            raise NotImplemented("no data-only output for cbar")
        if plot_type[-1] == "l":
            lvec = ["a", "c"]
        elif plot_type[-1] == "a":
            lvec = ["a"]
        else:
            lvec = ["c"]
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
                cbar_mat = data[pfx + 'cbar_sld_{l}'.format(l=l)][0]
                colors = cmap(cbar_mat.reshape(-1))
                collection[indx].set_color(colors)
                ttl.set_text('')
            out = [collection[i] for i in range(len(collection))]
            out.append(ttl)
            out = tuple(out)
            return out
        def animate(tind):
            for indx, l in enumerate(lvec):
                cbar_mat = data[pfx + 'cbar_sld_{l}'.format(l=l)][tind]
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
        if data_only:
            raise NotImplemented("no data-only output for bulkp")
        mpl.animation.Animation._blit_draw = _blit_draw
        fig, ax = plt.subplots()
        ax.set_xlabel('Position in electrode [{unit}]'.format(unit=Lunit))
        ax.set_ylabel('Potential of cathode [nondim]')
        ttl = ax.text(0.5, 1.05, ttl_fmt.format(perc=0),
                transform = ax.transAxes, verticalalignment="center",
                horizontalalignment="center")
        bulkp = pfx + 'phi_bulk_{l}'.format(l=l)
        Ltrode = dD['L'][l] * Lfac
        datay = data[bulkp][0]
        ymin = np.min(data[bulkp]) - 0.2
        ymax = np.max(data[bulkp]) + 0.2
        numy = len(datay)
        xmin = 0
        xmax = Ltrode
        datax = np.linspace(xmin, xmax, numy)
        ax.set_ylim((ymin, ymax))
        ax.set_xlim((xmin, xmax))
        # returns tuble of line objects, thus comma
        line1, = ax.plot(datax, datay)
        def init():
            line1.set_ydata(np.ma.array(datax, mask=True))
            ttl.set_text('')
            return line1, ttl
        def animate(tind):
            datay = data[bulkp][tind]
            line1.set_ydata(datay)
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
        ani.save("mpet_{type}.mp4".format(type=plot_type), fps=30)
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
