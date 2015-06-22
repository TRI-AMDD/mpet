import sys
import os
import datetime

import numpy as np
import scipy.io as sio
import matplotlib as mpl
mpl.use("TkAgg")
import matplotlib.pyplot as plt
import matplotlib.animation as manim
import matplotlib.collections as mcollect
import matplotlib.patches as mpatch

import mpet
import mpet_params_IO

def show_data(indir, plot_type, print_flag, save_flag, data_only):
    pfx = 'mpet.'
    ttl_fmt = "% = {perc:2.1f}"
    # Read in the simulation results and calcuations data
    dataFileName = "output_data.mat"
    dataFile = os.path.join(indir, dataFileName)
    data = sio.loadmat(dataFile)
    try:
        data[pfx + 'td'][0][0]
    except KeyError:
        pfx = ''
    # Read in the parameters used to define the simulation
    paramFileName = "input_params.cfg"
    paramFile = os.path.join(indir, paramFileName)
    IO = mpet_params_IO.mpetIO()
    P = IO.getConfig(paramFile)
    D = IO.getDictFromConfig(P)
    mod = mpet.modMPET("mpet", D=D)
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

    # Plot defaults
    axtickfsize = 18
    labelfsize = 20
    lwidth = 3.
    markersize = 10
    mpl.rcParams['xtick.labelsize'] = axtickfsize
    mpl.rcParams['ytick.labelsize'] = axtickfsize
    mpl.rcParams['font.size'] = labelfsize
    mpl.rcParams['legend.fontsize'] = labelfsize - 2
    mpl.rcParams['lines.linewidth'] = lwidth
    mpl.rcParams['lines.markersize'] = markersize

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
        print "Specified psd_mean [nm]:", np.array(D['psd_mean_ac'])[trodes]*1e9
        print "Specified psd_stddev [nm]:", np.array(D['psd_stddev_ac'])[trodes]*1e9
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

    def cbar(c_sld, l):
        Nij = len(c_sld)
        r_vec, volfrac_vec = mod.get_unit_solid_discr(
                D['solidShape_ac'][l],
                D['solidType_ac'][l], Nij)
        return np.sum(c_sld*volfrac_vec)
    def read_areas():
        fi = open('graphite_data/graphite_area_counts_cm2.txt', 'r')
        fi.readline() # heading
        data_t = []
        data_a1 = []
        data_a2 = []
        data_a3 = []
        for line in fi:
            line = line.strip().split()
            data_t.append(float(line[0]))
            data_a1.append(float(line[1]))
            data_a2.append(float(line[2]))
            data_a3.append(float(line[3]))
        return (data_t, data_a1, data_a2, data_a3)
    def get_images_times():
        namebase = 'circ_trans_cropped_G2_int-20130204-'
        imgfolder = os.path.join(os.getcwd(), 'graphite_data',
                'circ-trans-crop', 'colorMods')
        imgfilenames = os.listdir(imgfolder)
        namesTrueTimes = []
        for imgfilename in imgfilenames:
#            imgfilename = imgfilename.strip(namebase)
#            imgfilename = imgfilename.strip('.jpg')
            if imgfilename[:5] == "movie":
                continue
            imgtimestr = imgfilename[35:41]
            imgtime = datetime.datetime.strptime(imgtimestr, '%H%M%S')
            namesTrueTimes.append((imgfilename, imgtime))
        namesTrueTimes = sorted(namesTrueTimes, key=lambda nt : nt[1])
        t0 = namesTrueTimes[0][1]
        namesTimes = []
        for nTT in namesTrueTimes:
            imgfilename = nTT[0]
            imgtime = (nTT[1] - t0).total_seconds()
            namesTimes.append((imgfilename, imgtime))
        return namesTimes
    def get_image(tcur, namesTimes, imgdir):
        indx = 0
        for i, nT in enumerate(namesTimes):
            indx = i
            if tcur <= nT[1]:
                break
        imgname = namesTimes[indx][0]
        image = mpl.image.imread(os.path.join(imgdir, imgname))
        return image
    def smooth(vec, k=3):
        if k == 0:
            return vec
        if k % 2 == 0:
            raise Exception("smooth needs odd argument")
        ksides = (k-1)/2
        svec = np.zeros(len(vec))
        svec[ksides:-ksides] = [
                np.sum(vec[i-ksides:i+1+ksides])/float(k)
                for i in range(ksides, len(svec) - ksides)]
        svec[0] = vec[0]
        svec[-1] = vec[-1]
        for sind in range(1, ksides):
            ksides = sind
            svec[sind] = (np.sum(vec[:sind+1+ksides]) /
                    float(2*ksides+1))
            svec[-(sind+1)] = (np.sum(vec[-(sind+1)-ksides:]) /
                    float(2*ksides+1))
        return svec
    def flattenlist(nlist):
        flist = []
        for item in nlist:
            if type(item) in (list, tuple):
                flist.extend(flatten(item))
            else:
                flist.append(item)
        return flist

    # Plot voltage profile
    if plot_type == "v":
#        Vstd = 0.120
        voltage = Vstd - (k*Tref/e)*data[pfx + 'phi_applied'][0]
#        voltage = -data[pfx + 'phi_applied'][0]
        ffvec = data[pfx + 'ffrac_1'][0]
        from xlrd import open_workbook
        wb = open_workbook('graphite_data/ocp.xls')
        sheet = wb.sheet_by_index(0)
        exp_ff = sheet.col_values(0)[1:]
        exp_v = sheet.col_values(1)[1:]
        if data_only:
            return ffvec, voltage
        fig, ax = plt.subplots()
        ax.plot(ffvec, voltage,
                label="Model")
        ax.plot(exp_ff, exp_v, 'or', markersize=int(0.7*markersize),
                label="Expt. data")
#        ax.plot(times*td, voltage)
#        xmin = np.min(ffvec)
#        xmax = np.max(ffvec)
        xmin = 0.
        xmax = 1.
        ax.set_xlim((xmin, xmax))
        ax.legend(loc="best")
#        if not D['delPhiEqFit_ac'][1]:
##            ax.axhline(y=Vstd, xmin=xmin, xmax=xmax, linestyle='--', color='g')
#            ax.axhline(y=Vstd, linestyle='--', color='g')
        ax.set_xlabel("Cathode Filling Fraction [dimensionless]")
        ax.set_ylabel("Voltage [V]")
#        ax.set_ylim((Vstd_c - 0.3, Vstd_c + 0.4))
#        ax.set_ylim((2, 5))
        if save_flag:
            fig.savefig("mpet_v.png", bbox_inches="tight")
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
        scl = 0.5
        fig, ax = plt.subplots(Npart_ac[l], Nvol_ac[l], squeeze=False,
                sharey=True, figsize=(scl*16, scl*9))
#        lmbda = data[pfx + "lambda_ac"][0][l]
        k0 = D['k0_ac'][l]
        sol_c1_str_base = pfx + "c1_sld_trode{l}vol{{j}}part{{i}}".format(l=l)
        sol_c2_str_base = pfx + "c2_sld_trode{l}vol{{j}}part{{i}}".format(l=l)
        sol_p_str = pfx + "phi_{l}".format(l=l)
        lyte_c_str = pfx + "c_lyte_{l}".format(l=l)
        lyte_p_str = pfx + "phi_lyte_{l}".format(l=l)
        ylim = (0, 1.01)
        ffvec = data[pfx + 'ffrac_{l}'.format(l=l)][0]
        datax = times*td
#        datax = ffvec
        for i in range(Npart_ac[l]):
            for j in range(Nvol_ac[l]):
                sol1_str = sol_c1_str_base.format(i=i, j=j)
                sol2_str = sol_c2_str_base.format(i=i, j=j)
                # Remove axis ticks
#                ax[i, j].xaxis.set_major_locator(plt.NullLocator())
                csld1 = data[sol1_str]
                csld2 = data[sol2_str]
                c1bar_vec = [cbar(csld1[ii], l) for ii in range(len(times))]
                c2bar_vec = [cbar(csld2[ii], l) for ii in range(len(times))]
                c_lyte = data[lyte_c_str][:, j]
                act_O = c_lyte
                phi_lyte = data[lyte_p_str][:, j]
                voltage = Vstd - (k*Tref/e)*phi_lyte
                phi_m = data[sol_p_str][:, j]
                # PARAMETERS
                Omga = data[pfx + "Omga_{l}".format(l=l)][0][j][i]
                Omgb = data[pfx + 'Omgb_{l}'.format(l=l)][0][j][i]
                Omgc = data[pfx + 'Omgc_{l}'.format(l=l)][0][j][i]
                B = data[pfx + "B_ac"][0][l]
                kappa = data[pfx + "kappa_{l}".format(l=l)][0][j][i]
                EvdW = data[pfx + "EvdW_ac"][0][l]
                beta_s = data[pfx + "beta_s_{l}".format(l=l)][0][j][i]
                T = data[pfx + "T"][0][0]
                k0 = data[pfx + "k0_{l}".format(l=l)][0][j][i]
                alpha = data[pfx + "alpha_ac"][0][l]
                delta_L = data[pfx + "delta_L_{l}".format(l=l)][0][j][i]
                Ds = data[pfx + "Dsld_{l}".format(l=l)][0][j][i]
                Nij = len(csld1[0])
                r_vec, volfrac_vec = mod.get_unit_solid_discr(
                        D['solidShape_ac'][l],
                        D['solidType_ac'][l], Nij)
                dr = r_vec[1] - r_vec[0]
                mu1_R = np.zeros(csld1.shape)
                mu2_R = np.zeros(csld1.shape)
                for tind in range(len(times)):
                    (mu1_R[tind, :], mu2_R[tind, :]) = mod.calc_mu12_R(
                            csld1[tind], csld2[tind],
                            c1bar_vec[tind], c2bar_vec[tind],
                            Omga, Omgb, Omgc, B, kappa, EvdW,
                            beta_s, T)[0:2]
                eta1 = ((mu1_R[:, -1] + phi_m[:]) -
                        (T*np.log(act_O) + phi_lyte[:]))
                eta2 = ((mu2_R[:, -1] + phi_m[:]) -
                        (T*np.log(act_O) + phi_lyte[:]))
                gamma1_ts = (1./(csld1[:, -1]*(1-csld1[:, -1])))
                gamma2_ts = (1./(csld2[:, -1]*(1-csld2[:, -1])))
                ResistRxn = 1./(1./gamma1_ts + 1./gamma2_ts)
                ecd1 = (k0 * act_O**(1-alpha) *
                        np.exp(mu1_R[:,-1])**(1-alpha) / gamma1_ts)
                ecd2 = (k0 * act_O**(1-alpha) *
                        np.exp(mu2_R[:,-1])**(1-alpha) / gamma2_ts)
                rd1_vec = ecd1 * np.exp(-alpha*eta1/T)
                ox1_vec = ecd1 * np.exp((1-alpha)*eta1/T)
                rd2_vec = ecd2 * np.exp(-alpha*eta2/T)
                ox2_vec = ecd2 * np.exp((1-alpha)*eta2/T)
                rxn1_vec = mod.R_BV(k0, alpha, csld1[:, -1], csld2[:, -1],
                        Omgb, act_O, np.exp(mu1_R[:, -1]/T), eta1, T)
                rxn2_vec = mod.R_BV(k0, alpha, csld2[:, -1], csld1[:, -1],
                        Omgb, act_O, np.exp(mu2_R[:, -1]/T), eta2, T)
                rnet1_vec = rd1_vec - ox1_vec
                rnet2_vec = rd2_vec - ox2_vec
#                line1, = ax[i, j].plot(datax, -eta1, label=r"$-\eta_1$")
#                line2, = ax[i, j].plot(datax, -eta2, label=r"$-\eta_2$")
#                line1, = ax[i, j].plot(datax, mu1_R[:, -1],
#                        label=r"$\mu1_R$")
#                line2, = ax[i, j].plot(datax, mu2_R[:, -1],
#                        label=r"$\mu2_R$")
#                line1, = ax[i, j].plot(datax, phi_m, label=r"$\phi_m$")
#                line1, = ax[i, j].plot(datax, T*np.log(act_O), label=r"$\mu_O$")
#                line1, = ax[i, j].plot(datax, phi_lyte,
#                        label=r"$\phi_{lyte}$")
                line1, = ax[i, j].plot(datax, voltage, color="b",
                        label=r"$V_{appl}$")
                for tl in ax[i, j].get_yticklabels():
                    tl.set_color('b')
                ax[i, j].set_ylabel(r"Voltage [V]", color="b")
                ax[i, j].set_xlabel("Time [s]")
#                line1, = ax[i, j].plot(datax, 1./gamma1_ts,
#                        label=r"$1/\gamma_1$")
#                line2, = ax[i, j].plot(datax, 1./gamma2_ts,
#                        label=r"$1/\gamma_2$")
                ax2 = ax[i, j].twinx()
#                line1, = ax2.plot(datax, 1./gamma1_ts, linestyle='--',
#                        label=r"$1/\gamma_1$")
#                line2, = ax2.plot(datax, 1./gamma2_ts, linestyle='--',
#                        label=r"$1/\gamma_2$")
                line1, = ax2.plot(datax, ResistRxn, linestyle='--',
                        color='red',
                        label=r"$Resist_{tot}$")
                for tl in ax2.get_yticklabels():
                    tl.set_color('r')
#                ax2.set_ylabel("Reaction Resistance [nondimensional]",
                ax2.set_ylabel(r"$R_\mathrm{rxn}$ [nondimensional]",
                        color="r")
#                line1, = ax2.plot(datax, rxn1_vec, linestyle='--',
#                        label=r"$rxn_1$")
##                line3, = ax2.plot(datax, rd1_vec, linestyle=':',
##                        label=r"$rd_1$")
##                line4, = ax2.plot(datax, ox1_vec, linestyle='-.',
##                        label=r"$ox_1$")
#                line2, = ax2.plot(datax, rxn2_vec, linestyle='--',
#                        label=r"$rxn_2$")
##                line5, = ax2.plot(datax, rd2_vec, linestyle=':',
##                        label=r"$rd_2$")
##                line6, = ax2.plot(datax, ox2_vec, linestyle='-.',
##                        label=r"$ox_2$")
#                ax2.legend(loc="upper right")
#                ax[i, j].legend(loc="upper left")
                ax2.set_ylim((0., 100))
                ax[i, j].set_ylim((-0.05, 0.10))
##                # homog particles
##                mu_R = Omga*(1-2*csld) + 1*np.log(csld/(1-csld))
##                # ACR particles
##                sld_pt_indx = 50
##                Nij = csld.shape[1]
##                cbar = np.tile(np.sum(csld, axis=1)/Nij, [Nij, 1]).transpose()
##                cstmp = np.zeros((len(datax), Nij + 2))
##                cstmp[:, 1:-1] = csld
##                cstmp[0] = D['cwet_ac'][l]
##                cstmp[-1] = D['cwet_ac'][l]
##                dxs = 1./Nij
##                curv = np.diff(cstmp,2)/(dxs**2)
##                kappa = data[pfx + "kappa_{l}".format(l=l)][i, j]
##                B = data[pfx + "B_ac"][0][l]
##                mu_R = ( mu_reg_sln(csld, Omga) - kappa*curv
##                    + B*(csld - cbar) )
##                mu_R = mu_R[:, sld_pt_indx]
##                act_R = np.exp(mu_R)
##                eta = (mu_R + phi_m) - (mu_O + phi_lyte)
##                BVgamma_ts = 1./(1-csld[:, sld_pt_indx])
#                # diffn with delPhiEqFit
#                csurf = csld[:,-1]
#                import delta_phi_fits
#                dphi_eq_ref = data[pfx + "dphi_eq_ref_ac"][0][l]
#                fits = delta_phi_fits.DPhiFits(D)
#                phifunc = fits.materialData[D["material_ac"][l]]
#                dphi_eq = phifunc(csurf, dphi_eq_ref)
#                eta = phi_m - phi_lyte - dphi_eq
#                act_R = csurf
#                BVgamma_ts = 1./(1-csurf)
#                # General
#                act_O = c_lyte
#                mu_O = np.log(act_O)
##                eta = np.linspace(-20, 20, 70)
#                BValpha = D['alpha_ac'][l]
#                BVecd = ( k0 * act_O**(1-BValpha)
#                    * act_R**(BValpha) / BVgamma_ts )
#                BVrate = ( BVecd *
#                    (np.exp(-BValpha*eta/1) - np.exp((1-BValpha)*eta/1)) )
##                Malpha = 0.5*(1 + (1/lmbda) * np.log(c_lyte/csld))
##                Mecd = ( k0 *
##                    c_lyte**((3-2*Malpha)/4.) *
##                    csld**((1+2*Malpha)/4.) )
##                Meta2 = np.exp(-eta**2/(4.*1*lmbda))
##                Mrate = ( Mecd * np.exp(-eta**2/(4.*1*lmbda)) *
##                    (np.exp(-Malpha*eta/1) - np.exp((1-Malpha)*eta/1)) )
#                line, = ax[i, j].plot(datax, BVecd)
        if save_flag:
            fig.savefig("Rxn_out.pdf", bbox_inches="tight")
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
        mpl.animation.Animation._blit_draw = _blit_draw
        if data_only:
            raise NotImplemented("no data-only output for elytec/p")
        fig, ax = plt.subplots()
        if plot_type == "elytec":
            ymin = 0
            ymax = 2.2
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

    # Plot cslice
    elif plot_type in ["cslice_a", "cslice_c"]:
        fig, ax = plt.subplots()
        sol1 = data[pfx + "c1_sld_trode1vol0part0"]
        sol2 = data[pfx + "c2_sld_trode1vol0part0"]
        locinx_start = 0
        locinx_end = len(sol1[0, :]) - 1
        ax.plot(times*td, sol1[:, locinx_end])
#        ax.plot(times*td, sol2[:, locinx_end])
        print sol1[-1, locinx_start], sol2[-1, locinx_start]
        print sol1[-1, locinx_end], sol2[-1, locinx_end]
        return fig, ax

    # Plot all solid concentrations or potentials
    elif plot_type in ["csld_c", "csld_a", "phisld_a", "phisld_c",
            "csld_col_c", "csld_col_a"]:
        t0ind = len(times) - 1
#        t0ind = 300
        print "t =", times[t0ind]*td, "s"
        l = (0 if plot_type[-1] == "a" else 1)
        if data_only:
            raise NotImplemented("no data-only output for csld/phisld")
        scl = 0.5
        fig, ax = plt.subplots(Npart_ac[l], Nvol_ac[l], squeeze=False,
                sharey=True, figsize=(scl*16, scl*9))
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
#                ax[i, j].xaxis.set_major_locator(plt.NullLocator())
#                ax[i, j].yaxis.set_major_locator(plt.NullLocator())
#                print data[sol1[i, j]].shape
                datay1 = data[sol1[i, j]][t0ind]
                datay2 = data[sol2[i, j]][t0ind]
#                print data1.shape
                numy = len(datay1)
                datax = np.linspace(0, lens[i, j]*1e6, numy)
                ax[i, j].set_ylim(ylim)
                ax[i, j].set_xlim((0, lens[i, j]*1e6))
                ax[i, j].set_xlabel(r"$r$ [$\mu$m]")
                ax[i, j].set_ylabel(r"$c$")
                line1, = ax[i, j].plot(datax, datay1, label=r"$c_1$")
                line2, = ax[i, j].plot(datax, datay2, label=r"$c_2$")
                ax[i, j].legend(loc="best")
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
#                fig.savefig("csld_c.pdf", bbox_inches="tight")
#                plt.show()
#                zz

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
#        l = (0 if plot_type[-1] == "a" else 1)
        l = 1
        i = 0
        j = 0
        t0ind = 0
#        t0ind = 470
        t0ind = 330
        if data_only:
            raise NotImplemented("no data-only output for csld/phisld")
        # Define colors!
        to_red = 0.3
        to_yellow = 0.85
        rgb_rSd = 0.678
        rgb_gSd = 0.604
        rgb_bSd = 0.0
        rgb_rS2 = 0.970
        rgb_gS2 = 0.486
        rgb_bS2 = 0.0
        rgb_rS1 = 1.00
        rgb_gS1 = 0.973
        rgb_bS1 = 0.0
        cdict = {
                "red" : [(0.0, rgb_rSd, rgb_rSd),
                         (to_red, rgb_rSd, rgb_rS2),
                         (to_yellow, rgb_rS2, rgb_rS1),
                         (1.0, rgb_rS1, rgb_rS1)],
                "green" : [(0.0, rgb_gSd, rgb_gSd),
                           (to_red, rgb_gSd, rgb_gS2),
                           (to_yellow, rgb_gS2, rgb_gS1),
                           (1.0, rgb_gS1, rgb_gS1)],
                "blue" : [(0.0, rgb_bSd, rgb_bSd),
                          (to_red, rgb_bSd, rgb_bS2),
                          (to_yellow, rgb_bS2, rgb_bS1),
                          (1.0, rgb_bS1, rgb_bS1)]
                }
        cmap = mpl.colors.LinearSegmentedColormap(
                "discrete", cdict)

#        fig = plt.figure()
        figscale = 1.1
        fig = plt.figure(figsize=(figscale*16, figscale*9))
        axcirc = fig.add_axes([0.04, 0.04, 0.44, 0.40])
        axcirc.set_frame_on(False)
        axcirc.set_axis_off()
        axcirc.set_aspect('equal')
        axmovie = fig.add_axes([0.025, 0.51, 0.475, 0.470])
        axmovie.set_axis_off()
#        axmu1det = fig.add_axes([0.04, 0.04, 0.44, 0.40])
#        axmu2det = fig.add_axes([0.04, 0.56, 0.44, 0.40])
        axff = fig.add_axes([0.54, 0.56, 0.44, 0.40])
        axcsld = fig.add_axes([0.54, 0.08, 0.44, 0.40])
#        # tmp axes for cbar plot
#        fig = plt.figure(figsize=(figscale*8, figscale*9))
#        axff = fig.add_axes([0.04, 0.06, 0.88, 0.40])
#        axcsld = fig.add_axes([0.04, 0.56, 0.88, 0.40])
        csld1 = data[pfx + "c1_sld_trode1vol0part0"]
        csld2 = data[pfx + "c2_sld_trode1vol0part0"]
        phisld = data[pfx + "phi_1"]
        philyte = data[pfx + "phi_lyte_1"]
        datay1 = csld1[t0ind]
        datay2 = csld2[t0ind]
        # PARAMETERS
        Omga = data[pfx + "Omga_{l}".format(l=l)][0][j][i]
        Omgb = data[pfx + 'Omgb_{l}'.format(l=l)][0][j][i]
        Omgc = data[pfx + 'Omgc_{l}'.format(l=l)][0][j][i]
        B = data[pfx + "B_ac"][0][l]
        kappa = data[pfx + "kappa_{l}".format(l=l)][0][j][i]
        EvdW = data[pfx + "EvdW_ac"][0][l]
        beta_s = data[pfx + "beta_s_{l}".format(l=l)][0][j][i]
        T = data[pfx + "T"][0][0]
        k0 = data[pfx + "k0_{l}".format(l=l)][0][j][i]
        alpha = data[pfx + "alpha_ac"][0][l]
        delta_L = data[pfx + "delta_L_{l}".format(l=l)][0][j][i]
        Ds = data[pfx + "Dsld_{l}".format(l=l)][0][j][i]
        Nij = len(csld1[0])
        r_vec, volfrac_vec = mod.get_unit_solid_discr(
                D['solidShape_ac'][l],
                D['solidType_ac'][l], Nij)
        dr = r_vec[1] - r_vec[0]
        act_O = 1.
        c1bar_vec = [cbar(csld1[i], l) for i in range(len(times))]
        c2bar_vec = [cbar(csld2[i], l) for i in range(len(times))]
        mu1_R = np.zeros(csld1.shape)
        mu2_R = np.zeros(csld1.shape)
        mu1_reg = np.zeros(csld1.shape)
        mu2_reg = np.zeros(csld1.shape)
        mu1_B = np.zeros(csld1.shape)
        mu2_B = np.zeros(csld1.shape)
        mu1_vdW = np.zeros(csld1.shape)
        mu2_vdW = np.zeros(csld1.shape)
        mu1_intr = np.zeros(csld1.shape)
        mu2_intr = np.zeros(csld1.shape)
        Flux1 = np.zeros((csld1.shape[0], csld1.shape[1]+1))
        Flux2 = np.zeros((csld1.shape[0], csld1.shape[1]+1))
        # Chemical Potential
        for tind in range(len(times)):
            (mu1_R[tind, :], mu2_R[tind, :],
                    mu1_reg[tind], mu2_reg[tind],
                    mu1_B[tind], mu2_B[tind],
                    mu1_vdW[tind], mu2_vdW[tind],
                    mu1_intr[tind], mu2_intr[tind]) = mod.calc_mu12_R(
                    csld1[tind], csld2[tind],
                    c1bar_vec[tind], c2bar_vec[tind],
                    Omga, Omgb, Omgc, B, kappa, EvdW, beta_s, T)
        mu1_kappa = mu1_R - (mu1_reg + mu1_B + mu1_vdW + mu1_intr)
        mu2_kappa = mu2_R - (mu2_reg + mu2_B + mu2_vdW + mu2_intr)
        # Overpotential
        eta1 = ((mu1_R[:, -1] + phisld[:, -1]) -
                (T*np.log(act_O) + philyte[:, -1]))
        eta2 = ((mu2_R[:, -1] + phisld[:, -1]) -
                (T*np.log(act_O) + philyte[:, -1]))
        # Reactions
        rxn1_vec = mod.R_BV(k0, alpha, csld1[:, -1], csld2[:, -1],
                Omgb, act_O, np.exp(mu1_R[:, -1]/T), eta1, T)
        rxn2_vec = mod.R_BV(k0, alpha, csld2[:, -1], csld1[:, -1],
                Omgb, act_O, np.exp(mu2_R[:, -1]/T), eta2, T)
        Flux1_bc_vec = 0.5*delta_L*rxn1_vec
        Flux2_bc_vec = 0.5*delta_L*rxn2_vec
        for tind in range(len(times)):
            Flux1[tind, :], Flux2[tind, :] = mod.calc_Flux12(
                    csld1[tind], csld2[tind],
                    mu1_R[tind], mu2_R[tind], Ds,
                    Flux1_bc_vec[tind], Flux2_bc_vec[tind],
                    dr, T)
        rxn1flux_vec = Flux1[:, -1]
        rxn2flux_vec = Flux2[:, -1]
        dataybar = 0.5*(datay1 + datay2)
        smcount = 7
        numy = len(datay1)
        p_len = psd_len_ac[l][0][0, 0]
        datax = np.linspace(0, p_len, numy)
        ylim = (0, 1)
#        t_offset = d_t[0]
        t_offset = 0

        # csld plot
        axcsld.set_ylim(ylim)
        axcsld.set_xlim((0, p_len*1e6))
        line1, = axcsld.plot(datax*1e6, datay1,
                color="blue",
                label=r'$c_1$')
        line2, = axcsld.plot(datax*1e6, datay2,
                color="green",
                label=r'$c_2$')
        line3, = axcsld.plot(datax*1e6, 0.5*(datay1 + datay2),
                color='red', linestyle='--',
                label=r'$\overline{c} = 0.5*(c_1 + c_2)$')
        axcsld.legend(loc='best')
        axcsld.set_xlabel(r"$r$ [$\mu$m]")
        axcsld.set_ylabel(r"$c$")

#        # csld plot (axff axes)
#        axff.set_ylim(ylim)
#        axff.set_xlim((0, p_len*1e6))
#        line1, = axff.plot(datax*1e6, datay1,
#                color="blue",
#                label=r'$c_1$')
#        line2, = axff.plot(datax*1e6, datay2,
#                color="green",
#                label=r'$c_2$')
#        line3, = axff.plot(datax*1e6, 0.5*(datay1 + datay2),
#                color='red', linestyle='--',
#                label=r'$\overline{c} = 0.5*(c_1 + c_2)$')
#        axff.legend(loc='best')
#        axff.set_xlabel(r"$r$ [$\mu$m]")
#        axff.set_ylabel(r"$c$")
#        bar1indx = 120
#        bar2indx = 273
#        axff.fill_between(datax[0:bar1indx+1]*1e6, 0, 1,
#                color=cmap(0.1), alpha=0.7)
#        axff.fill_between(datax[bar1indx:bar2indx+1]*1e6, 0, 1,
#                color=cmap(0.6), alpha=0.7)
#        axff.fill_between(datax[bar2indx:]*1e6, 0, 1,
#                color=cmap(0.9), alpha=0.7)

#        # csld plot schematic-like
#        axcsld.set_ylim((0, 1))
#        axcsld.set_xlabel(r"$r$ [$\mu$m]")
#        axcsld.yaxis.set_ticks([])
#        axcsld.set_axis_bgcolor('none')
#        axspines = axcsld.spines
#        for spine in ['top', 'bottom', 'left', 'right']:
#            axspines[spine].set_visible(False)
#        axcsld.xaxis.set_ticks_position('bottom')
#        xmax = p_len*1e6
#        axcsld.set_xlim((0, xmax))
#        # Draw graphite planes
#        gpatches = []
#        ytop = 0.8
#        gthick = 0.015
#        ng = 9
#        if not ng % 2: raise Exception("ng must be odd")
#        gylocs = np.linspace(0.02, ytop-gthick/2., ng)
#        for yloc in gylocs:
#            gpatches.append(mpatch.Rectangle((0, yloc-gthick/2.),
#                xmax, gthick, facecolor='black', edgecolor='black',
#                ))
#        gcollect = mcollect.PatchCollection(gpatches,
#                match_original=True)
#        axcsld.add_collection(gcollect)
#        # Color discretization
#        datar = datax*1e6
#        dr = datar[1] - datar[0]
#        # Set up bounding pie slice
#        b0 = 0.3
#        a0 = p_len*1e6
#        nx = 100
#        xpietop = 18./24. * a0
#        xclip_r = np.linspace(xpietop, xmax, nx)
#        yclip_r = b0*np.sqrt(1 - (xclip_r/a0)**2) + ytop
#        # Colored ellipses on top
#        # Set up a clipping path for the ellipses
#        path_array = np.array([
#            [xmax, ytop],
#            [0, ytop],
#            ])
#        xy_r = np.hstack((xclip_r.reshape(-1,1),
#                          yclip_r.reshape(-1,1)))
#        path_array = np.vstack((path_array, xy_r))
#        clippath = mpl.path.Path(path_array)
#        clippatch = mpatch.PathPatch(clippath, facecolor='none',
#                edgecolor='none')
#        axcsld.add_patch(clippatch)
#        # Add ellipses which are the "top-down" color on top of the
#        # pie slice. Have to add these from outer to inner so they can
#        # stack on each other...
#        ellipfills = []
#        ellipfcol = cmap(dataybar[0])
#        ellipfills.append(mpatch.Ellipse((0, ytop),
#            2*xmax, 2*xmax*b0/a0,
#            facecolor=ellipfcol, linewidth=0.))
#        nelip = numy
#        for indxtop in range(1, nelip):
#            r_outer = xmax - dr/2. - (indxtop-1)*dr
#            ellipfills.append(mpatch.Ellipse((0, ytop),
#                2*r_outer, 2*r_outer*b0/a0,
#                facecolor=ellipfcol, linewidth=0))
#        ellipcollect = mcollect.PatchCollection(ellipfills,
#                match_original=True)
#        axcsld.add_collection(ellipcollect)
#        ellipcollect.set_clip_path(clippatch)
#        # Colored rectangles inside
#        nlayers = ng - 1
#        nA = nB = nlayers/2
#        nrect = numy
#        rectsA = []
#        rectsB = []
#        for inx in range(nA):
#            rectsA.append([])
#            rectsB.append([])
##        cmapcsld = plt.get_cmap('gray')
##        cmapcsld = plt.get_cmap('Greys')
##        cmapcsld = plt.get_cmap('cool')
#        cmapcsld = plt.get_cmap('summer')
#        col1 = cmapcsld(datay1[0])
#        col2 = cmapcsld(datay2[0])
#        alph1 = 1.0
#        alph2 = 1.0
#        h = gylocs[-1] - gylocs[-2] - gthick
#        ysA = []
#        ysB = []
#        for inx in range(nA):
#            ysA.append(gylocs[-(2 + 2*inx)] + gthick/2.)
#            ysB.append(gylocs[-(3 + 2*inx)] + gthick/2.)
#        for inx in range(nA):
#            yA = ysA[inx]
#            yB = ysB[inx]
#            rectsA[inx].append(mpatch.Rectangle((0, yA), dr/2., h,
#            facecolor=col1, edgecolor=col1, alpha=alph1,
#            linewidth=0,
#            ))
#            rectsB[inx].append(mpatch.Rectangle((0, yB), dr/2., h,
#            facecolor=col1, edgecolor=col1, alpha=alph1,
#            linewidth=0,
#            ))
#        for indxrect in range(1, nrect-1):
#            ri = datar[indxrect] - dr/2.
#            for indxlyr in range(nA):
#                yA = ysA[indxlyr]
#                yB = ysB[indxlyr]
#                rectsA[indxlyr].append(mpatch.Rectangle((ri, yA), dr, h,
#                facecolor=col1, edgecolor=col1, alpha=alph1,
#                linewidth=0,
#                ))
#                rectsB[indxlyr].append(mpatch.Rectangle((ri, yB), dr, h,
#                facecolor=col1, edgecolor=col1, alpha=alph1,
#                linewidth=0,
#                ))
#        ri = datar[-1] - dr/2.
#        for indxlyr in range(nA):
#            yA = ysA[indxlyr]
#            yB = ysB[indxlyr]
#            rectsA[indxlyr].append(mpatch.Rectangle((ri, yA), dr/2., h,
#            facecolor=col1, edgecolor=col1, alpha=alph1,
#            linewidth=0,
#            ))
#            rectsB[indxlyr].append(mpatch.Rectangle((ri, yB), dr/2., h,
#            facecolor=col1, edgecolor=col1, alpha=alph1,
#            linewidth=0,
#            ))
#        rectcollsA = []
#        rectcollsB = []
#        for indxlyr in range(nA):
#            rectcollsA.append(mcollect.PatchCollection(rectsA[indxlyr],
#                match_original=True,
#                ))
#            rectcollsB.append(mcollect.PatchCollection(rectsB[indxlyr],
#                match_original=True,
#                ))
#            axcsld.add_collection(rectcollsA[indxlyr])
#            axcsld.add_collection(rectcollsB[indxlyr])
#        csldcolors1 = cmapcsld(datay1)
#        csldcolors2 = cmapcsld(datay2)
#        for indxlyr in range(nA):
#            rectcollsA[indxlyr].set_color(csldcolors1)
#            rectcollsB[indxlyr].set_color(csldcolors2)
#        topcol = cmap(smooth(dataybar, smcount))
#        ellipcollect.set_color(topcol[::-1])

        # colored circle plot
        axcirc.set_ylim(ylim)
        datar = datax/datax[-1]/2.05 # normalized to 1/2
        dr = datar[1] - datar[0]
        ncirc = numy
        circfills = np.empty(ncirc, dtype=object)
        circfills = []
        col = cmap(dataybar[0])
        circfills.append(mpatch.Wedge((0.5, 0.5), dr/2., 0, 360,
            width=dr/2.,
            facecolor=col, edgecolor=col,
            ))
        for indxcirc in range(1, ncirc-1):
            ri = datar[indxcirc] - dr/2.
            ro = datar[indxcirc] + dr/2.
            circfills.append(mpatch.Wedge((0.5, 0.5), ro, 0, 360,
                width=dr,
                facecolor=col, edgecolor=col,
                ))
        circfills.append(mpatch.Wedge((0.5, 0.5), datar[-1], 0, 360,
            width=dr/2.,
            facecolor=col, edgecolor=col,
            ))
        collection = mcollect.PatchCollection(circfills,
                match_original=True,
                )
        axcirc.add_collection(collection)
        colors = cmap(smooth(dataybar, smcount))
        collection.set_color(colors)

#        # Experimental images movie
#        namesTimes = get_images_times()
#        imgdir = os.path.join(os.getcwd(), 'graphite_data',
#                'circ-trans-crop', 'colorMods')
#        image = get_image(times[t0ind]*td, namesTimes, imgdir)
#        img = axmovie.imshow(image)
#        patch = mpatch.Circle((178, 162), radius=135, transform=axmovie.transData)
#        img.set_clip_path(patch)

#        # ff or soc plot
#        ffvec = data[pfx + 'ffrac_{l}'.format(l=l)][0]
#        axff.plot(times*td, ffvec)
#        axff.set_xlabel("time [s]")
#        axff.set_ylabel("Filling Fraction")
#        ffcirc, = axff.plot(times[0]*td, ffvec[0], 'or',
#                markerfacecolor='none',
#                markeredgecolor='red',
#                markeredgewidth=2,
#                )

        # voltage plot on exp images axis
        axmovie.set_axis_on()
        axmovie.set_position([0.08, 0.56, 0.40, 0.40])
        ffvec = data[pfx + 'ffrac_{l}'.format(l=l)][0]
        tvec = times*td
        voltage = Vstd - (k*Tref/e)*data[pfx + 'phi_applied'][0]
        axmovie.plot(ffvec, voltage)
#        axmovie.plot(tvec, voltage)
        axmovie.set_xlabel("Cathode Filling Fraction [dimensionless]")
        axmovie.set_ylabel("Voltage [V]")
        ffcirc, = axmovie.plot(ffvec[t0ind], voltage[t0ind], 'or',
#        ffcirc, = axmovie.plot(tvec[t0ind], voltage[t0ind], 'or',
                markerfacecolor='none',
                markeredgecolor='red',
                markeredgewidth=2,
                )

#        # areas plot
#        d_t, d_a1, d_a2, d_a3 = read_areas()
#        sc = 1e0
#        axff.plot(d_t, sc*np.array(d_a1), marker='o',
#                linestyle='None',
#                markerfacecolor=(rgb_rS1, rgb_gS1, rgb_bS1),
#                markeredgecolor=(rgb_rS1, rgb_gS1, rgb_bS1),
#                label="expt: Stage 1")
#        axff.plot(d_t, sc*np.array(d_a2), marker='s',
#                linestyle='None',
#                markerfacecolor=(rgb_rS2, rgb_gS2, rgb_bS2),
#                markeredgecolor=(rgb_rS2, rgb_gS2, rgb_bS2),
#                label="expt: Stage 2")
#        axff.plot(d_t, sc*np.array(d_a3), marker='^',
#                linestyle='None',
#                markerfacecolor=(rgb_rSd, rgb_gSd, rgb_bSd),
#                markeredgecolor=(rgb_rSd, rgb_gSd, rgb_bSd),
#                label="expt: Dilute")
#        d_Atot = np.array(d_a1) + np.array(d_a2) + np.array(d_a3)
#        axff.set_xlabel(r'Time (s)')
#        axff.set_ylabel(r'Area (cm$^2$)')
#        axff.yaxis.major.formatter.set_powerlimits((0,0)) # sci. not'n
#        area_calcs = np.zeros((len(times), 3))
#        for tind in range(len(times)):
#            cbar = 0.5*(csld1[tind, :] + csld2[tind, :])
##            cbar = smooth(0.5*(csld1[tind, :] + csld2[tind, :]),
##                    smcount)
#            s3ind = np.where(cbar < to_red)
#            s2ind = np.where(
#                    np.logical_and(cbar > to_red, cbar < to_yellow))
#            s1ind = np.where(cbar > to_yellow)
#            Atot = np.pi*D['psd_mean_ac'][l]**2 * 1e2**2 # cm^2
#            # note for cylinder, volfrac = areafrac
#            area_calcs[tind, 0] = Atot*np.sum(volfrac_vec[s1ind])
#            area_calcs[tind, 1] = Atot*np.sum(volfrac_vec[s2ind])
#            area_calcs[tind, 2] = Atot*np.sum(volfrac_vec[s3ind])
#        axff.plot(times*td + t_offset, sc*area_calcs[:, 0],
#                linestyle='-',
#                color=(rgb_rS1, rgb_gS1, rgb_bS1),
#                label="sim: Stage 1")
#        axff.plot(times*td + t_offset, sc*area_calcs[:, 1],
#                linestyle='-',
#                color=(rgb_rS2, rgb_gS2, rgb_bS2),
#                label="sim: Stage 2")
#        axff.plot(times*td + t_offset, sc*area_calcs[:, 2],
#                linestyle='-',
#                color=(rgb_rSd, rgb_gSd, rgb_bSd),
#                label="sim: Dilute")
##        axff.plot(d_t, d_Atot, '.k')
##        axff.plot(times*td, np.sum(area_calcs, axis=1), '-k')
#        ffline = axff.axvline(times[t0ind]*td + t_offset,
#                linestyle='--',
#                linewidth=0.7*lwidth,
#                color='#808080')
#        axff.legend(loc='best')

#        # ff --> rxns plot
#        rxn1vec = data[pfx + 'rxn1'][0]
#        rxn2vec = data[pfx + 'rxn2'][0]
#        axff.plot(times*td, rxn1vec)
#        axff.plot(times*td, rxn2vec)
#        axff.set_xlabel("time [s]")
#        axff.set_ylabel("Layer Rxn Rate [dimless]")
#        ffline = axff.axvline(times[0]*td + t_offset)
#        # ff --> rxns plot from flux calcs
#        axff.plot(times*td, rxn1flux_vec)
#        axff.plot(times*td, rxn2flux_vec)
#        axff.set_xlabel("time [s]")
#        axff.set_ylabel("Layer Rxn Rate [dimless]")
#        axff.set_ylim((0, np.max(rxn1flux_vec.max(),
#            rxn2flux_vec.max())))
#        ffline = axff.axvline(times[0]*td)
#        # Flux plot
#        edges = p_len*1e6*np.hstack((0, (r_vec[0:-1] + r_vec[1:])/2., 1))
#        urdata1 = Flux1
#        urdata2 = Flux2
#        urline1, = axff.plot(edges, urdata1[t0ind, :])
#        urline2, = axff.plot(edges, urdata2[t0ind, :])
#        axff.set_xlim((0, p_len*1e6))
##        axff.set_ylim((min(np.nanmin(Flux1), np.nanmin(Flux2)),
##            max(np.nanmax(Flux1), np.nanmaxFlux2))))
#        axff.set_ylim((-0.012, 0.012))
#        axff.set_ylabel(r"Flux")
        # mu_R plot
        edges = p_len*1e6*np.hstack((0, (r_vec[0:-1] + r_vec[1:])/2., 1))
        urdata1 = mu1_R
        urdata2 = mu2_R
        urline1, = axff.plot(r_vec*p_len*1e6, mu1_R[t0ind, :])
        urline2, = axff.plot(r_vec*p_len*1e6, mu2_R[t0ind, :])
        axff.set_xlim((0, p_len*1e6))
        axff.set_ylim((min(np.nanmin(mu1_R), np.nanmin(mu2_R)),
            max(np.nanmax(mu1_R), np.nanmax(mu2_R))))
        axff.set_ylabel(r"$\mu_R$")
#        # dc/dt plot
#        edges = p_len*1e6*np.hstack((0, (r_vec[0:-1] + r_vec[1:])/2., 1))
#        dcdt1 = np.diff(Flux1*2*np.pi*edges)/volfrac_vec
#        dcdt2 = np.diff(Flux2*2*np.pi*edges)/volfrac_vec
#        urdata1 = dcdt1
#        urdata2 = dcdt2
#        urline1, = axff.plot(r_vec*p_len*1e6, urdata1[0, :])
#        urline2, = axff.plot(r_vec*p_len*1e6, urdata2[0, :])
#        axff.set_xlim((0, p_len*1e6))
##        axff.set_ylim((min(np.nanmin(dcdt1), np.nanmin(dcdt2)),
##            max(np.nanmax(dcdt1), np.nanmax(dcdt2))))
#        axff.set_ylim((-10, 10))
#        axff.set_ylabel(r"dcdt")

#        # D_eff Calculations
##        # for t0ind = 490
##        rdil_0 = 0
##        rdil_1 = 7.6
##        rstg2_0 = 8.9
##        rstg2_1 = 18.6
##        rstg1_0 = 19.5
##        rstg1_1 = datax[-1]*1e6
##        t210 = rstg2_0
##        t211 = 11.0
##        t220 = 11.8
##        t221 = 17.1
##        t230 = 18.11
##        t231 = rstg2_1
##        # t0ind = 290
##        rdil_0 = 0
##        rdil_1 = 15
##        rstg2_0 = 16.6
##        rstg2_1 = 21.4
##        rstg1_0 = 23.2
##        rstg1_1 = datax[-1]*1e6
##        t210 = rstg2_0
##        t211 = 11.0
##        t220 = 11.8
##        t221 = 17.1
##        t230 = 18.11
##        t231 = rstg2_1
#        # t0ind = 590
#        rdil_0 = 0
#        rdil_1 = 3.2
#        rstg2_0 = 4.3
#        rstg2_1 = 16.0
#        rstg1_0 = 17.3
#        rstg1_1 = datax[-1]*1e6
#        t210 = rstg2_0
#        t211 = 5.4
#        t220 = 6.3
#        t221 = 10.9
#        t230 = 11.9
#        t231 = rstg2_1
#        print "Time:", times[t0ind]*td, "s"
#        rdil_ind0 = np.argmin(np.abs(datax*1e6 - rdil_0))
#        rdil_ind1 = np.argmin(np.abs(datax*1e6 - rdil_1))
#        rstg2_ind0 = np.argmin(np.abs(datax*1e6 - rstg2_0))
#        rstg2_ind1 = np.argmin(np.abs(datax*1e6 - rstg2_1))
#        rstg1_ind0 = np.argmin(np.abs(datax*1e6 - rstg1_0))
#        rstg1_ind1 = np.argmin(np.abs(datax*1e6 - rstg1_1))
#        ti210 = np.argmin(np.abs(datax*1e6 - t210))
#        ti211 = np.argmin(np.abs(datax*1e6 - t211))
#        ti220 = np.argmin(np.abs(datax*1e6 - t220))
#        ti221 = np.argmin(np.abs(datax*1e6 - t221))
#        ti230 = np.argmin(np.abs(datax*1e6 - t230))
#        ti231 = np.argmin(np.abs(datax*1e6 - t231))
##        axcsld.axvline(datax[ti210]*1e6, color='magenta')
##        axcsld.axvline(datax[ti211]*1e6, color='magenta')
##        axcsld.axvline(datax[ti220]*1e6, color='magenta')
##        axcsld.axvline(datax[ti221]*1e6, color='magenta')
##        axcsld.axvline(datax[ti230]*1e6, color='magenta')
##        axcsld.axvline(datax[ti231]*1e6, color='magenta')
#        axff.axvline(datax[rdil_ind0]*1e6, color='k')
#        axff.axvline(datax[rdil_ind1]*1e6, color='k')
#        axff.axvline(datax[rstg2_ind0]*1e6, color='r')
#        axff.axvline(datax[rstg2_ind1]*1e6, color='r')
#        axff.axvline(datax[rstg1_ind0]*1e6, color='y')
#        axff.axvline(datax[rstg1_ind1]*1e6, color='y')
#        axcsld.axvline(datax[rdil_ind1]*1e6, color='k')
#        axcsld.axvline(datax[rdil_ind1]*1e6, color='k')
#        axcsld.axvline(datax[rstg2_ind0]*1e6, color='r')
#        axcsld.axvline(datax[rstg2_ind1]*1e6, color='r')
#        axcsld.axvline(datax[rstg1_ind0]*1e6, color='y')
#        axcsld.axvline(datax[rstg1_ind1]*1e6, color='y')
#        # Dilute:
#        print "Dilute"
#        F1avg_dil = np.mean(Flux1[rdil_ind0:rdil_ind1])
#        F2avg_dil = np.mean(Flux2[rdil_ind0:rdil_ind1])
#        gradCavg_dil = ((dataybar[rdil_ind1] - dataybar[rdil_ind0]) /
#                (r_vec[rdil_ind1] - r_vec[rdil_ind0]))
#        Dndim_dil = (F1avg_dil + F2avg_dil)/gradCavg_dil
#        D_dil = Dndim_dil*p_len**2/td
#        print D_dil
#        # Stage 2
#        print "Stage 2"
#        F1avg_stg2 = np.mean(Flux1[rstg2_ind0:rstg2_ind1])
#        F2avg_stg2 = np.mean(Flux2[rstg2_ind0:rstg2_ind1])
#        gradCavg_stg2 = ((dataybar[rstg2_ind1] - dataybar[rstg2_ind0]) /
#                (r_vec[rstg2_ind1] - r_vec[rstg2_ind0]))
#        Dndim_stg2 = (F1avg_stg2 + F2avg_stg2)/gradCavg_stg2
#        D_stg2 = Dndim_stg2*p_len**2/td
#        print D_stg2
#        # Stage 1
#        print "Stage 1"
#        F1avg_stg1 = np.mean(Flux1[rstg1_ind0:rstg1_ind1])
#        F2avg_stg1 = np.mean(Flux2[rstg1_ind0:rstg1_ind1])
#        gradCavg_stg1 = ((dataybar[rstg1_ind1] - dataybar[rstg1_ind0]) /
#                (r_vec[rstg1_ind1] - r_vec[rstg1_ind0]))
#        Dndim_stg1 = (F1avg_stg1 + F2avg_stg1)/gradCavg_stg1
#        D_stg1 = Dndim_stg1*p_len**2/td
#        print D_stg1
#        # Stage 2 intermmediates
#        print "Stage 2, section 1"
#        F1t21 = np.mean(Flux1[ti210:ti211])
#        F2t21 = np.mean(Flux2[ti210:ti211])
#        gradC21 = ((dataybar[ti211] - dataybar[ti210]) /
#                (r_vec[ti211] - r_vec[ti210]))
#        print (F1t21 + F2t21)/gradC21 * p_len**2/td
#        print "Stage 2, section 2"
#        F1t22 = np.mean(Flux1[ti220:ti221])
#        F2t22 = np.mean(Flux2[ti220:ti221])
#        gradC22 = ((dataybar[ti221] - dataybar[ti220]) /
#                (r_vec[ti221] - r_vec[ti220]))
#        print (F1t22 + F2t22)/gradC22 * p_len**2/td
#        print "Stage 2, section 3"
#        F1t23 = np.mean(Flux1[ti230:ti231])
#        F2t23 = np.mean(Flux2[ti230:ti231])
#        gradC23 = ((dataybar[ti231] - dataybar[ti230]) /
#                (r_vec[ti231] - r_vec[ti230]))
#        print (F1t23 + F2t23)/gradC23 * p_len**2/td

#        # mu details
#        lldata1, lldata2, lldata3, lldata4, lldata5 = (
#                mu1_reg, mu1_B, mu1_vdW, mu1_kappa, mu1_intr)
#        lldata = [lldata1, lldata2, lldata3, lldata4, lldata5]
#        llline1, = axmu1det.plot(r_vec*p_len*1e6, lldata1[0, :], 'b--')
#        llline2, = axmu1det.plot(r_vec*p_len*1e6, lldata2[0, :], 'b-')
#        llline3, = axmu1det.plot(r_vec*p_len*1e6, lldata3[0, :], 'bx')
#        llline4, = axmu1det.plot(r_vec*p_len*1e6, lldata4[0, :], 'b:')
#        llline5, = axmu1det.plot(r_vec*p_len*1e6, lldata5[0, :], 'b*')
#        lllines = [llline1, llline2, llline3, llline4, llline5]
#        axmu1det.set_xlim((0, p_len*1e6))
#        axmu1det.set_ylim((-3, 3))
#        #
#        uldata1, uldata2, uldata3, uldata4, uldata5 = (
#                mu2_reg, mu2_B, mu2_vdW, mu2_kappa, mu2_intr)
#        uldata = [uldata1, uldata2, uldata3, uldata4, uldata5]
#        ulline1, = axmu2det.plot(r_vec*p_len*1e6, uldata1[0, :], 'g--')
#        ulline2, = axmu2det.plot(r_vec*p_len*1e6, uldata2[0, :], 'g-')
#        ulline3, = axmu2det.plot(r_vec*p_len*1e6, uldata3[0, :], 'gx')
#        ulline4, = axmu2det.plot(r_vec*p_len*1e6, uldata4[0, :], 'g:')
#        ulline5, = axmu2det.plot(r_vec*p_len*1e6, uldata5[0, :], 'g*')
#        ullines = [ulline1, ulline2, ulline3, ulline4, ulline5]
#        axmu2det.set_xlim((0, p_len*1e6))
#        axmu2det.set_ylim((-3, 3))

#        # Flux ul
#        edges = p_len*1e6*np.hstack((0, (r_vec[0:-1] + r_vec[1:])/2., 1))
#        uldata1 = Flux1
#        uldata2 = Flux2
#        ulline1, = axmu2det.plot(edges, uldata1[0, :])
#        ulline2, = axmu2det.plot(edges, uldata2[0, :])
#        axmu2det.set_xlim((0, p_len*1e6))
##        axff.set_ylim((min(np.nanmin(Flux1), np.nanmin(Flux2)),
##            max(np.nanmax(Flux1), np.nanmaxFlux2))))
#        axmu2det.set_ylim((-0.012, 0.012))
#        axmu2det.set_ylabel(r"Flux")

#        fig.savefig('tmp.png', bbox_inches='tight')
#        plt.show()

        def init():
            toblit = []
            datay1 = csld1[0]
            datay2 = csld2[0]
            dataybar = 0.5*(datay1 + datay2)
            dataybar = smooth(0.5*(datay1 + datay2), smcount)
            # csld
            line1.set_ydata(np.ma.array(csld1[0], mask=True))
            line2.set_ydata(np.ma.array(csld2[0], mask=True))
            toblit.extend([line1, line2])
            line3.set_ydata(np.ma.array(csld2[0], mask=True))
            toblit.extend([line3])
#            # csld schematic-like
#            csldcolors1 = cmapcsld(datay1)
#            csldcolors2 = cmapcsld(datay2)
#            for indxlyr in range(nA):
#                rectcollsA[indxlyr].set_color(csldcolors1)
#                rectcollsB[indxlyr].set_color(csldcolors2)
#                toblit.extend([rectcollsA[indxlyr],
#                               rectcollsB[indxlyr]])
#            topcol = cmap(dataybar)
#            ellipcollect.set_color(topcol[::-1])
#            toblit.append(ellipcollect)
#            # ff
#            ffcirc.set_xdata(np.ma.array(0, mask=True))
#            ffcirc.set_ydata(np.ma.array(0, mask=True))
#            toblit.append(ffcirc)
            # voltage
            ffcirc.set_xdata(np.ma.array(0, mask=True))
            ffcirc.set_ydata(np.ma.array(0, mask=True))
            toblit.append(ffcirc)
#            # areas plot
#            ffline.set_xdata(np.ma.array([0, 0], mask=True))
#            toblit.append(ffline)
            # flux, mu_R, dcdt
            ny = len(urline1.get_ydata())
            urline1.set_ydata(np.ma.array([0]*ny, mask=True))
            urline2.set_ydata(np.ma.array([0]*ny, mask=True))
            toblit.extend([urline1, urline2])
            # color circle
            colors = cmap(dataybar)
            collection.set_color(colors)
            toblit.append(collection)
#            # experimental images
#            toblit.append(img)
#            # mu1 details
#            for line in lllines:
#                ny = len(line.get_ydata())
#                line.set_ydata(np.ma.array([0]*ny, mask=True))
#                toblit.append(line)
#            # mu2 details
#            for line in ullines:
#                ny = len(line.get_ydata())
#                line.set_ydata(np.ma.array([0]*ny, mask=True))
#                toblit.append(line)
            # Flux ul
#            ny = len(ulline1.get_ydata())
#            ulline1.set_ydata(np.ma.array([0]*ny, mask=True))
#            ulline2.set_ydata(np.ma.array([0]*ny, mask=True))
#            toblit.extend(ulline1, ulline2)
            return tuple(toblit)

        def animate(tind):
#            tind += 500
            toblit = []
            datay1 = csld1[tind]
            datay2 = csld2[tind]
            dataybar_raw = 0.5*(datay1 + datay2)
            dataybar = smooth(0.5*(datay1 + datay2), smcount)
            # csld
            line1.set_ydata(datay1)
            line2.set_ydata(datay2)
            toblit.extend([line1, line2])
            line3.set_ydata(dataybar_raw)
            toblit.extend([line3])
#            # csld schematic-like
#            csldcolors1 = cmapcsld(datay1)
#            csldcolors2 = cmapcsld(datay2)
#            for indxlyr in range(nA):
#                rectcollsA[indxlyr].set_color(csldcolors1)
#                rectcollsB[indxlyr].set_color(csldcolors2)
#                toblit.extend([rectcollsA[indxlyr],
#                               rectcollsB[indxlyr]])
#            topcol = cmap(dataybar)[::-1]
#            ellipcollect.set_color(topcol)
#            toblit.append(ellipcollect)
#            # ff
#            ffcirc.set_xdata(times[tind]*td)
#            ffcirc.set_ydata(ffvec[tind])
#            toblit.append(ffcirc)
            # voltage
            ffcirc.set_xdata(ffvec[tind])
#            ffcirc.set_xdata(tvec[tind])
            ffcirc.set_ydata(voltage[tind])
            toblit.append(ffcirc)
#            # areas plot
#            ffline.set_xdata(2*[times[tind]*td + t_offset])
#            toblit.append(ffline)
            # flux, mu_R, dcdt
            urline1.set_ydata(urdata1[tind, :])
            urline2.set_ydata(urdata2[tind, :])
            toblit.extend([urline1, urline2])
            # color circle
            colors = cmap(dataybar)
            collection.set_color(colors)
            toblit.append(collection)
#            # experimental images movie
#            tcur = times[tind]*td + t_offset
#            image = get_image(tcur, namesTimes, imgdir)
#            img.set_array(image)
#            toblit.append(img)
#            # mu1 details
#            for indx in range(len(lllines)):
#                line = lllines[indx]
#                data = lldata[indx]
#                line.set_ydata(data[tind, :])
#                toblit.append(line)
#            # mu2 details
#            for indx in range(len(ullines)):
#                line = ullines[indx]
#                data = uldata[indx]
#                line.set_ydata(data[tind, :])
#                toblit.append(line)
#            # Flux ul
#            ulline1.set_ydata(uldata1[tind, :])
#            ulline2.set_ydata(uldata2[tind, :])
#            toblit.extend([ulline1, ulline2])
            return tuple(toblit)

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
                cbar_mat = 0.5*(
                        data[pfx + 'c1bar_sld_{l}'.format(l=l)][0] +
                        data[pfx + 'c2bar_sld_{l}'.format(l=l)][0])
                colors = cmap(cbar_mat.reshape(-1))
                collection[indx].set_color(colors)
                ttl.set_text('')
            out = [collection[i] for i in range(len(collection))]
            out.append(ttl)
            out = tuple(out)
            return out
        def animate(tind):
            for indx, l in enumerate(lvec):
                cbar_mat = 0.5*(
                        data[pfx + 'c1bar_sld_{l}'.format(l=l)][tind] +
                        data[pfx + 'c2bar_sld_{l}'.format(l=l)][tind])
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
        mpl.animation.Animation._blit_draw = _blit_draw
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
        ani.save("mpet_{type}.mp4".format(type=plot_type),
                fps=25, bitrate=5500,
#                writer='alzkes',
#                savefig_kwargs={'bbox_inches' : 'tight'},
                )
#                extra_args=['-vcodec', 'libx264'])

    return ani

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
