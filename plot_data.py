import sys
import os
import time

import numpy as np
import scipy.io as sio
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import matplotlib.animation as manim
import matplotlib.collections as mcollect
import matplotlib.patches as mpatch

import mpet_params_IO

def show_data(indir, plot_type, print_flag, save_flag, data_only):
    pfx = 'mpet.'
    ttl_fmt = "% = {perc:2.1f}"
    # Read in the simulation results and calcuations data
    dataFileName = "output_data.mat"
    dataFile = os.path.join(indir, dataFileName)
    data = sio.loadmat(dataFile)
    # Read in the parameters used to define the simulation
    paramFileName = "output_params.cfg"
    paramFile = os.path.join(indir, paramFileName)
    IO = mpet_params_IO.mpetIO()
    P = IO.getConfig(paramFile)
    D = IO.getDictFromConfig(P)
    # simulated (porous) electrodes
    Nvol_ac = np.array([D['Nvol_ac'][0], D['Nvol_ac'][1]])
    if Nvol_ac[0] >= 1: # If porous anode
        trodes = [0, 1]
    else:
        trodes = [1]
#    Ntrode = data[pfx + 'NumTrode']
    Ntrode = 2
    # Pick out some useful constants/calculated values
    k = D['k']                      # Boltzmann constant, J/(K Li)
    Tref = D['Tref']                # Temp, K
    e = D['e']                      # Charge of proton, C
    td = data[pfx + 'td'][0][0]     # diffusive time
    Vstd_ac = np.zeros(Ntrode)
    Nvol_ac = np.zeros(Ntrode, dtype=np.integer)
    Npart_ac = np.zeros(Ntrode, dtype=np.integer)
    solidType_ac = np.empty(Ntrode, dtype=object)
    solidShape_ac = np.empty(Ntrode, dtype=object)
    rxnType_ac = np.empty(Ntrode, dtype=object)
    psd_len_ac = np.empty(Ntrode, dtype=object)
    for l in trodes:
        Vstd_ac[l] = D['Vstd_ac'][l]
        # Replace the standard potential if a fit voltage curve was used.
        # Use the value that the simulation used in initialization.
        if D['delPhiEqFit_ac'][l]:
            Vstd_ac[l] = data[pfx + 'dphi_eq_ref_ac'][0][l]*(k*Tref/e)
        Nvol_ac[l] = D['Nvol_ac'][l]
        Npart_ac[l] = D['Npart_ac'][l]
        # Simulation type
        solidType_ac[l] = D['solidType_ac'][l]
        solidShape_ac[l] = D['solidShape_ac'][l]
        rxnType_ac[l] = D['rxnType_ac'][l]
        psd_len_ac[l] = data[pfx + "psd_lengths_{l}".format(l=l)]
    Vstd = Vstd_ac[1] - Vstd_ac[0]
    Nvol_s = D['Nvol_s']
    # Discretization (and associated porosity)
    dxc = data[pfx + "L_ac"][0][1]/Nvol_ac[1]
    dxvec = np.array(Nvol_ac[1] * [dxc])
    porosvec = np.array(Nvol_ac[1] * [data[pfx + "poros_ac"][0][1]])
    L_c = D['L_ac'][1] * 1e6
    cellsvec = dxc*np.arange(Nvol_ac[1]) + dxc/2.
    if Nvol_s:
        dxs = data[pfx + "L_s"][0][0]/Nvol_s
        dxvec_s = np.array(Nvol_s * [dxs])
        dxvec = np.hstack((dxvec_s, dxvec))
        poros_s = np.array(Nvol_s * [data[pfx + "poros_s"][0][0]])
        porosvec = np.hstack((poros_s, porosvec))
        L_s = D['L_s'] * 1e6
        cellsvec += L_s/L_c
        cellsvec_s = dxs*np.arange(Nvol_s) + dxs/2.
        cellsvec = np.hstack((cellsvec_s, cellsvec))
    if 0 in trodes:
        dxa = data[pfx + "L_ac"][0][0]/Nvol_ac[0]
        dxvec_a = np.array(Nvol_ac[0] * [dxa])
        dxvec = np.hstack((dxvec_a, dxvec))
        poros_a = np.array(Nvol_ac[0] * [data[pfx + "poros_ac"][0][0]])
        porosvec = np.hstack((poros_a, porosvec))
        L_a = D['L_ac'][0] * 1e6
        cellsvec += L_a/L_c
        cellsvec_a = dxa*np.arange(Nvol_ac[0]) + dxa/2.
        cellsvec = np.hstack((cellsvec_a, cellsvec))
    cellsvec *= L_c
    # Extract the reported simulation times
    times = data[pfx + 'phi_applied_times'][0]
    numtimes = len(times)
    tmin = np.min(times)
    tmax = np.max(times)
    # Simulation type
    profileType = D['profileType']
    # Colors for plotting concentrations
    to_yellow = 0.3
    to_red = 0.7

    # Print relevant simulation info
    if print_flag:
#        for i in trodes:
#            print "PSD_{l}:".format(l=l)
#            print psd_len_ac[l][0].transpose()
#            print "Actual psd_mean [nm]:", np.mean(psd_len_ac[l])
#            print "Actual psd_stddev [nm]:", np.std(psd_len_ac[l])
        print "profileType:", profileType
        print "Cell structure:"
        print (("porous anode | " if Nvol_ac[0] else "flat anode | ")
                + ("sep | " if Nvol_s else "") + "porous cathode")
        if Nvol_ac[0]:
            print "capacity ratio cathode:anode, 'z':", data[pfx + 'z'][0][0]
        print "solidType:", solidType_ac[trodes]
        print "solidShape", solidShape_ac[trodes]
        print "rxnType:", rxnType_ac[trodes]
        if profileType == "CC":
            print "C_rate:", D['Crate']
        else: # CV
            print "Vset:", D['Vset']
        print "Specified psd_mean [nm]:", np.array(D['mean_ac'])[trodes]*1e9
        print "Specified psd_stddev [nm]:", np.array(D['stddev_ac'])[trodes]*1e9
#        print "reg sln params:"
#        print data[pfx + "a"][0]
        if Nvol_s: print "Nvol_s:", Nvol_s
        print "Nvol_c:", Nvol_ac[1]
        if Nvol_ac[0]: print "Nvol_a:", Nvol_ac[0]
        print "Npart_c:", Npart_ac[1]
        if Nvol_ac[0]: print "Npart_a:", Npart_ac[0]
        print "Dp [m^2/s]:", D['Dp']
        print "Dm [m^2/s]:", D['Dm']
        print "Damb [m^2/s]:", data[pfx + "dim_Damb"][0][0]
        print "td [s]:", data[pfx + "td"][0][0]
        print "k0 [A/m^2]:", np.array(D['k0_ac'])[trodes]
        for l in trodes:
            trode = ("a" if l == 0 else "c")
            if rxnType_ac[l] == "BV":
                print "alpha_" + trode + ":", D['alpha_ac'][l]
            elif rxnType_ac[l] == "Marcus":
                print "lambda_" + trode + "/(kTref):", data[pfx + "lmbda_ac"][0][l]
            elif rxnType_ac[l] == "MHC":
                print "MHC_Aa_" + trode + ":", data[pfx +
                        "MHC_Aa_{l}".format(l=l)][0][0]
                print "b_" + trode + ":", data[pfx + "MHC_b_ac"][0][l]
            if D['simBulkCond_ac'][l]:
                print (trode + " bulk conductivity loss: Yes -- " +
                        "dim_mcond [S/m]: " + str(D['mcond_ac'][l]))
            else:
                print trode + " bulk conductivity loss: No"
            if D['simSurfCond_ac'][l]:
                print (trode + " surface conductivity loss: Yes -- " +
                        "dim_scond [S]: " + str(D['scond_ac'][l]))
            else:
                print trode + " surface conductivity loss: No"

    # Plot voltage profile
    if plot_type == "v":
        voltage = Vstd - (k*Tref/e)*data[pfx + 'phi_applied'][0]
#        voltage = -data[pfx + 'phi_applied'][0]
        ffvec = data[pfx + 'ffrac_1'][0]
        if data_only:
            return ffvec, voltage
        fig, ax = plt.subplots()
        ax.plot(ffvec, voltage)
#        xmin = np.min(ffvec)
#        xmax = np.max(ffvec)
        xmin = 0.
        xmax = 1.
        ax.set_xlim((xmin, xmax))
#        if not D['delPhiEqFit']:
##            ax.axhline(y=Vstd_c, xmin=xmin, xmax=xmax, linestyle='--', color='g')
#            ax.axhline(y=Vstd_c, linestyle='--', color='g')
        ax.set_xlabel("Cathode Filling Fraction [dimensionless]")
        ax.set_ylabel("Voltage [V]")
#        ax.set_ylim((Vstd_c - 0.3, Vstd_c + 0.4))
#        ax.set_ylim((2, 5))
        if save_flag:
            fig.savefig("mpet_v.png")
        return fig, ax

    # Plot surface conc.
    if plot_type in ["surf_c", "surf_a"]:
        l = (0 if plot_type[-1] == "a" else 1)
        if data_only:
            raise NotImplemented("no data-only output for surf")
        fig, ax = plt.subplots(Npart_ac[l], Nvol_ac[l], squeeze=False,
                sharey=True)
        str_base = pfx + "c_sld_trode{l}vol{{j}}part{{i}}".format(l=l)
        ylim = (0, 1.01)
        datax = times
        for i in range(Npart_ac[l]):
            for j in range(Nvol_ac[l]):
                sol_str = str_base.format(i=i, j=j)
                # Remove axis ticks
                ax[i, j].xaxis.set_major_locator(plt.NullLocator())
                datay = data[sol_str][:,-1]
                line, = ax[i, j].plot(times, datay)
        return fig, ax

    # Plot misc stuff about reactions
    if plot_type in ["rxnp_c", "rxnp_a"]:
        l = (0 if plot_type[-1] == "a" else 1)
        if data_only:
            raise NotImplemented("no data-only output for rxnp")
        fig, ax = plt.subplots(Npart_ac[l], Nvol_ac[l], squeeze=False,
                sharey=True)
#        lmbda = data[pfx + "lambda_ac"][0][l]
        k0 = D['k0_ac'][l]
        sol_c_str_base = pfx + "c_sld_trode{l}vol{{j}}part{{i}}".format(l=l)
        sol_p_str = pfx + "phi_{l}".format(l=l)
        lyte_c_str = pfx + "c_lyte_{l}".format(l=l)
        lyte_p_str = pfx + "phi_lyte_{l}".format(l=l)
        ylim = (0, 1.01)
        ffvec = data[pfx + 'ffrac_{l}'.format(l=l)][0]
        datax = times
#        datax = ffvec
        for i in range(Npart_ac[l]):
            for j in range(Nvol_ac[l]):
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
                BValpha = D['alpha_ac'][l]
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
            fig.savefig("Rxn_out.png")
        return fig, ax

    # Plot SoC profile
    if plot_type in ["soc_c", "soc_a"]:
        l = (0 if plot_type[-1] == "a" else 1)
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
            fig.savefig("mpet_soc.png")
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
        anode = pfx + 'c_lyte_0'
        cath = pfx + 'c_lyte_1'
        ax.set_xlabel('Battery Position [um]')
        cvec = data[cath]
        if Nvol_s:
            cvec_s = data[sep]
            cvec = np.hstack((cvec_s, cvec))
        if 0 in trodes:
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
        ffvec = data[pfx + 'ffrac_1'][0]
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
            fig.savefig("mpet_current.png")
        return fig, ax

    # Plot electrolyte concentration or potential
    elif plot_type == "elytec" or plot_type == "elytep":
        matplotlib.animation.Animation._blit_draw = _blit_draw
        if data_only:
            raise NotImplemented("no data-only output for elytec/p")
        fig, ax = plt.subplots()
        if plot_type == "elytec":
            ymin = 0
            ymax = 1.5
            ax.set_ylabel('Concentration of electrolyte [nondim]')
            sep = pfx + 'c_lyte_s'
            anode = pfx + 'c_lyte_0'
            cath = pfx + 'c_lyte_1'
        elif plot_type == "elytep":
            ymin = -50
            ymax = 50
            ax.set_ylabel('Potential of electrolyte [nondim]')
            sep = pfx + 'phi_lyte_s'
            anode = pfx + 'phi_lyte_0'
            cath = pfx + 'phi_lyte_1'
        ax.set_xlabel('Battery Position [um]')
        ttl = ax.text(0.5, 1.05, ttl_fmt.format(perc=0),
                transform = ax.transAxes, verticalalignment="center",
                horizontalalignment="center")
        datay = data[cath][0]
        L_c = D['L_ac'][1] * 1e6
        Ltot = L_c
        if Nvol_s:
            datay_s = data[sep][0]
            datay = np.hstack((datay_s, datay))
            L_s = D['L_s'] * 1e6
            Ltot += L_s
        else:
            L_s = 0
        if 0 in trodes:
            datay_a = data[anode][0]
            datay = np.hstack((datay_a, datay))
            L_a = D['L_ac'][0] * 1e6
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
#        line1, = ax.plot(datax, datay)
        line1, = ax.plot(datax, datay, '-')
        ax.axvline(x=L_a, linestyle='--', color='g')
        ax.axvline(x=(L_a+L_s), linestyle='--', color='g')
        def init():
            line1.set_ydata(np.ma.array(datax, mask=True))
            ttl.set_text('')
            return line1, ttl
        def animate(tind):
            datay = data[cath][tind]
            if Nvol_s:
                datay_s = data[sep][tind]
                datay = np.hstack((datay_s, datay))
            if 0 in trodes:
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
        l = (0 if plot_type[-1] == "a" else 1)
        if data_only:
            raise NotImplemented("no data-only output for csld/phisld")
        fig, ax = plt.subplots(Npart_ac[l], Nvol_ac[l], squeeze=False,
                sharey=True)
        sol1 = np.empty((Npart_ac[l], Nvol_ac[l]), dtype=object)
        sol2 = np.empty((Npart_ac[l], Nvol_ac[l]), dtype=object)
        lens = np.zeros((Npart_ac[l], Nvol_ac[l]))
        lines = np.empty((Npart_ac[l], Nvol_ac[l], 2), dtype=object)
        fills = np.empty((Npart_ac[l], Nvol_ac[l], 3), dtype=object)
        if plot_type in ["csld_a", "csld_col_a"]:
            str_base = pfx + "c_sld_trode0vol{j}part{i}"
            ylim = (0, 1.01)
        elif plot_type in ["csld_c", "csld_col_c"]:
            str_base1 = pfx + "c1_sld_trode1vol{j}part{i}"
            str_base2 = pfx + "c2_sld_trode1vol{j}part{i}"
            ylim = (0, 1.01)
        elif plot_type == "phisld_a": # plot_type == "phisld"
            str_base = pfx + "p_sld_trode0vol{j}part{i}"
            ylim = (-10, 20)
        elif plot_type == "phisld_c":
            str_base = pfx + "p_sld_trode1vol{j}part{i}"
            ylim = (-10, 20)
        for i in range(Npart_ac[l]):
            for j in range(Nvol_ac[l]):
                sol1[i, j] = str_base1.format(i=i, j=j)
                sol2[i, j] = str_base2.format(i=i, j=j)
                lens[i, j] = psd_len_ac[l][0][j, i]
                # Remove axis ticks
                ax[i, j].xaxis.set_major_locator(plt.NullLocator())
#                ax[i, j].yaxis.set_major_locator(plt.NullLocator())
                datay1 = data[sol1[i, j]][0]
                datay2 = data[sol2[i, j]][0]
                numy = len(datay1)
                datax = np.linspace(0, lens[i, j], numy)
                ax[i, j].set_ylim(ylim)
                ax[i, j].set_xlim((0, lens[i, j]))
                line1, = ax[i, j].plot(datax, datay1)
                line2, = ax[i, j].plot(datax, datay2)
                lines[i, j, 0] = line1
                lines[i, j, 1] = line2
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
            for i in range(Npart_ac[l]):
                for j in range(Nvol_ac[l]):
                    datax = np.zeros(data[sol1[i, j]][0].shape)
                    lines[i, j, 0].set_ydata(np.ma.array(datax, mask=True))
                    lines[i, j, 1].set_ydata(np.ma.array(datax, mask=True))
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
            for i in range(Npart_ac[l]):
                for j in range(Nvol_ac[l]):
                    datay1 = data[sol1[i, j]][tind]
                    datay2 = data[sol2[i, j]][tind]
                    lines[i, j, 0].set_ydata(datay1)
                    lines[i, j, 1].set_ydata(datay2)
                    datax = lines[i, j, 0].get_xdata()
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

    # Plot all solid concentrations or potentials
    elif plot_type in ["csldcirc"]:
        l = (0 if plot_type[-1] == "a" else 1)
        if data_only:
            raise NotImplemented("no data-only output for csld/phisld")
        # Define colors!
        to_red = 0.3
        to_yellow = 0.6
        cdict = {
                "red" : [(0.0, 0.0, 0.0),
                         (to_red, 0.0, 0.7),
                         (to_yellow, 0.7, 0.7),
                         (1.0, 0.7, 0.7)],
                "green" : [(0.0, 0.0, 0.0),
                           (to_red, 0.0, 0.0),
                           (to_yellow, 0.0, 0.6),
                           (1.0, 0.6, 0.6)],
                "blue" : [(0.0, 0.4, 0.4),
                          (to_red, 0.0, 0.0),
                          (to_yellow, 0.0, 0.0),
                          (1.0, 0.0, 0.0)]
                }
        cmap = matplotlib.colors.LinearSegmentedColormap(
                "discrete", cdict)
        fig = plt.figure()
        axcirc = fig.add_axes([0.04, 0.04, 0.44, 0.94])
        axcirc.set_frame_on(False)
        axcirc.set_axis_off()
        axcirc.set_aspect('equal')
        axff = fig.add_axes([0.54, 0.56, 0.44, 0.40])
        axcsld = fig.add_axes([0.54, 0.08, 0.44, 0.40])
        sol1 = data[pfx + "c1_sld_trode1vol0part0"]
        sol2 = data[pfx + "c2_sld_trode1vol0part0"]
        datay1 = sol1[0]
        datay2 = sol2[0]
        dataybar = 0.5*(datay1 + datay2)
        numy = len(datay1)
        p_len = psd_len_ac[l][0][0, 0]
        datax = np.linspace(0, p_len, numy)
        ylim = (0, 1)
        axcirc.set_ylim(ylim)
        axcsld.set_ylim(ylim)
#        axcirc.set_rlim((0, p_len))
        axcsld.set_xlim((0, p_len*1e9))
        # csld plot
        line1, = axcsld.plot(datax*1e9, datay1)
        line2, = axcsld.plot(datax*1e9, datay2)
        axcsld.set_xlabel("r [nm]")
        axcsld.set_ylabel(r"$\widetilde{c}$")
        datar = datax/datax[-1]/2.05 # normalized to 1/2
        dr = datar[1] - datar[0]
        ncirc = numy
        circfills = np.empty(ncirc, dtype=object)
        circfills = []
        col = cmap(dataybar[0])
        circfills.append(mpatch.Wedge((0.5, 0.5), dr/2., 0, 360,
            width=dr/2.,
            facecolor=col, edgecolor=col,
#            transform=axcirc.transAxes,
            ))
        for indxcirc in range(1, ncirc-1):
            ri = datar[indxcirc] - dr/2.
            ro = datar[indxcirc] + dr/2.
            circfills.append(mpatch.Wedge((0.5, 0.5), ro, 0, 360,
                width=dr,
                facecolor=col, edgecolor=col,
#                transform=axcirc.transAxes,
                ))
        circfills.append(mpatch.Wedge((0.5, 0.5), datar[-1], 0, 360,
            width=dr/2.,
            facecolor=col, edgecolor=col,
#            transform=axcirc.transAxes,
            ))
        collection = mcollect.PatchCollection(circfills,
                match_original=True,
                )
        axcirc.add_collection(collection)
        # ff or soc plot
        ffvec = data[pfx + 'ffrac_{l}'.format(l=l)][0]
        axff.plot(times*td, ffvec)
        axff.set_xlabel("time [s]")
        axff.set_ylabel("Filling Fraction")
        ffcirc, = axff.plot(times[0]*td, ffvec[0], 'or',
                markersize=7, markerfacecolor='none',
                markeredgecolor='red',
                markeredgewidth=2,
                )
#        # Ghetto-rigged v2
#        circfills[0] = axcirc.fill_between(
#                np.linspace(0.0, 2*np.pi, 100),
#                0*np.ones(100), dr/2.*np.ones(100),
#                facecolor=col, edgecolor=col,
#                )
#        circfills[-1] = axcirc.fill_between(
#                np.linspace(0.0, 2*np.pi, 100),
#                datar[-1]*np.ones(100),
#                (datar[-1] - dr/2.)*np.ones(100),
#                facecolor=col, edgecolor=col,
#                )
#        for indxcirc in range(1, ncirc-1):
#            ri = datar[indxcirc] - dr/2.
#            ro = datar[indxcirc] + dr/2.
#            circfills[indxcirc] = axcirc.fill_between(
#                    np.linspace(0.0, 2*np.pi, 100),
#                    ri*np.ones(100), ro*np.ones(100),
#                    facecolor=col, edgecolor=col,
#                    )
#        # Ghetto-rigged circle v1...
#        wherevec = np.linspace(0, 1, 100)
#        wherevec = np.sin(datar/p_len/(0.1))
#        nslice = 100
#        nthet = nslice + 1
#        thetavec = np.linspace(0, 2*np.pi, nthet)
#        slicefills = np.empty((nslice, 3), dtype=object)
#        for indxthet in range(nslice):
#            thet1 = thetavec[indxthet]
#            thet2 = thetavec[indxthet+1]
#            slicefills[indxthet, 0] = axcirc.fill_betweenx(
#                datar, thet1, thet2,
#                facecolor='yellow', edgecolor='none',
#                where=(wherevec > 0.8),
#                )
#            slicefills[indxthet, 1] = axcirc.fill_betweenx(
#                datar, thet1, thet2,
#                facecolor='red', edgecolor='none',
#                where=((wherevec > 0.3) & (wherevec < 0.8)),
#                )
#            slicefills[indxthet, 2] = axcirc.fill_betweenx(
#                datar, thet1, thet2,
#                facecolor='blue', edgecolor='none',
#                where=(wherevec < 0.3),
#                )
        def init():
#            datax = np.zeros(data[sol1[i, j]][0].shape)
            datay1 = sol1[0]
            datay2 = sol2[0]
            dataybar = 0.5*(datay1 + datay2)
            line1.set_ydata(np.ma.array(datay1, mask=True))
            line2.set_ydata(np.ma.array(datay2, mask=True))
            ffcirc.set_xdata(np.ma.array(0, mask=True))
            ffcirc.set_ydata(np.ma.array(0, mask=True))
            colors = cmap(dataybar)
#            print type(collection)
#            print dir(collection)
            collection.set_color(colors)
            return tuple([collection, line1, line2, ffcirc])
        def animate(tind):
#            datax = np.zeros(data[sol1[i, j]][tind].shape)
            datay1 = sol1[tind]
            datay2 = sol2[tind]
            dataybar = 0.5*(datay1 + datay2)
            line1.set_ydata(datay1)
            line2.set_ydata(datay2)
            ffcirc.set_xdata(times[tind]*td)
            ffcirc.set_ydata(ffvec[tind])
            colors = cmap(dataybar)
            collection.set_color(colors)
            return tuple([collection, line1, line2, ffcirc])

    # Plot average solid concentrations
    elif plot_type in ["cbar_c", "cbar_a", "cbar_full"]:
        if data_only:
            raise NotImplemented("no data-only output for cbar")
        if plot_type[-1] == "l":
            lvec = [0, 1]
        elif plot_type[-1] == "a":
            lvec = [0]
        else:
            lvec = [1]
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
            cmap = matplotlib.colors.LinearSegmentedColormap(
                    "discrete", cdict)
        # Smooth colormap changes:
        if color_changes == "smooth":
#            cmap = matplotlib.cm.RdYlGn_r # A default green-yellow-red map
            # generated with colormap.org
            cmaps = np.load("colormaps_custom.npz")
            cmap_data = cmaps["GnYlRd_3"]
            cmap = matplotlib.colors.ListedColormap(cmap_data/255.)

        # Implement hack to be able to animate title
        matplotlib.animation.Animation._blit_draw = _blit_draw
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
            psd_len = psd_len_ac[l][0]
            len_max = np.max(psd_len)
            len_min = np.min(psd_len)
            ax.patch.set_facecolor('white')
            # Don't stretch axes to fit figure -- keep 1:1 x:y ratio.
            ax.set_aspect('equal', 'box')
            # Don't show axis ticks
            ax.xaxis.set_major_locator(plt.NullLocator())
            ax.yaxis.set_major_locator(plt.NullLocator())
            ax.set_xlim(0, 1.)
            ax.set_ylim(0, float(Npart_ac[l])/Nvol_ac[l])
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
            spacing = 1.0 / Nvol_ac[l]
            size_fracs = 0.4*np.ones((Nvol_ac[l], Npart_ac[l]))
            if len_max != len_min:
                size_fracs = (psd_len - len_min)/(len_max - len_min)
            sizes = (size_fracs * (1 - size_frac_min) + size_frac_min) / Nvol_ac[l]
            # Create rectangle "patches" to add to figure axes.
            rects = np.empty((Nvol_ac[l], Npart_ac[l]), dtype=object)
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
        l = (0 if plot_type[-1] == "a" else 1)
        if data_only:
            raise NotImplemented("no data-only output for bulkp")
        matplotlib.animation.Animation._blit_draw = _blit_draw
        fig, ax = plt.subplots()
        ymin = -1
        ymax = 10
        ax.set_xlabel('Position in electrode [um]')
        ax.set_ylabel('Potential of cathode [nondim]')
        ttl = ax.text(0.5, 1.05, ttl_fmt.format(perc=0),
                transform = ax.transAxes, verticalalignment="center",
                horizontalalignment="center")
        bulkp = pfx + 'phi_{l}'.format(l=l)
        Ltrode = D['L_ac'][l] * 1e6
        datay = data[bulkp][0]
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
        ani.save("mpet_{type}.mp4".format(type=plot_type), fps=20)
#                extra_args=['-vcodec', 'libx264'])

    return ani

# This is a block of code which messes with some matplotlib internals
# to allow for animation of a title. See
# http://stackoverflow.com/questions/17558096/animated-title-in-matplotlib
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
