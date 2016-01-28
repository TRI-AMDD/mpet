import os
import ConfigParser
import pickle

import numpy as np

import delta_phi_fits
import elyte_CST

class mpetIO():

    def getConfigs(self, paramfile="params.cfg"):
        # system-level config
        P_s = ConfigParser.RawConfigParser()
        P_s.optionxform = str
        P_s.read(paramfile)

        # electrode config(s)
        P_e = {}
        # cathode config
        paramfile_c = P_s.get('Electrodes', 'cathode')
        # Look in the same directory as the input paramfile, wherever
        # that is.
        paramfileLoc = os.path.split(paramfile)[0]
        paramfile_c = os.path.join(paramfileLoc, paramfile_c)
        P_c = ConfigParser.RawConfigParser()
        P_c.optionxform = str
        P_c.read(paramfile_c)
        P_e["c"] = P_c

        # anode config
        if P_s.getint('Sim Params', 'Nvol_a') >= 1:
            paramfile_a = P_s.get('Electrodes', 'anode')
            paramfile_a = os.path.join(paramfileLoc, paramfile_a)
            P_a = ConfigParser.RawConfigParser()
            P_a.optionxform = str
            P_a.read(paramfile_a)
            P_e["a"] = P_a
        return P_s, P_e

    def getDictsFromConfigs(self, P_s, P_e):
        # Dictionary of dimensional system and electrode parameters
        dD_s = {}
        dD_e = {}
        # Dictionary of nondimensional system and electrode parameters
        ndD_s = {}
        ndD_e = {}

        # General particle classification
        ndD_s["2varTypes"] = ["diffn2", "CHR2", "homog2", "homog2_sdn"]
        ndD_s["1varTypes"] = ["ACR", "diffn", "CHR", "homog",
                "homog_sdn"]

        # Simulation Parameters
        ndD_s["profileType"] = P_s.get('Sim Params', 'profileType')
        dD_s["Crate"] = P_s.getfloat('Sim Params', 'Crate')
        try:
            dD_s["Vmax"] = P_s.getfloat('Sim Params', 'Vmax')
        except ConfigParser.NoOptionError:
            dD_s["Vmax"] = 1e10
        try:
            dD_s["Vmin"] = P_s.getfloat('Sim Params', 'Vmin')
        except ConfigParser.NoOptionError:
            dD_s["Vmin"] = -1e10
        dD_s["Vset"] = P_s.getfloat('Sim Params', 'Vset')
        ndD_s["capFrac"] = P_s.getfloat('Sim Params', 'capFrac')
        dD_s["tend"] = P_s.getfloat('Sim Params', 'tend')
        ndD_s["tsteps"] = P_s.getfloat('Sim Params', 'tsteps')
        Tabs = dD_s["Tabs"] = P_s.getfloat('Sim Params', 'T')
        Rser = dD_s["Rser"] = P_s.getfloat('Sim Params', 'Rser')
        ndD_s["Nvol"] = {"a": P_s.getint('Sim Params', 'Nvol_a'),
                       "c": P_s.getint('Sim Params', 'Nvol_c'),
                       "s": P_s.getint('Sim Params', 'Nvol_s')}
        ndD_s["Npart"] = {"a": P_s.getint('Sim Params', 'Npart_a'),
                        "c": P_s.getint('Sim Params', 'Npart_c')}
        ndD_s["trodes"] = ["c"]
        if ndD_s["Nvol"]["a"] >= 1:
            ndD_s["trodes"].append("a")
        dD_s["k0_foil"] = P_s.getfloat('Electrodes', 'k0_foil')

        # Particle info
        dD_s["psd_mean"] = {"a": P_s.getfloat('Particles', 'mean_a'),
                          "c":  P_s.getfloat('Particles', 'mean_c')}
        dD_s["psd_stddev"] = {"a": P_s.getfloat('Particles', 'stddev_a'),
                            "c": P_s.getfloat('Particles', 'stddev_c')}
        ndD_s["cs0"] = {"a": P_s.getfloat('Particles', 'cs0_a'),
                "c": P_s.getfloat('Particles', 'cs0_c')}

        # Conductivity
        ndD_s["simBulkCond"] = {"a": P_s.getboolean('Conductivity', 'simBulkCond_a'),
                "c": P_s.getboolean('Conductivity', 'simBulkCond_c')}
        dD_s["mcond"] = {"a": P_s.getfloat('Conductivity', 'mcond_a'),
                "c": P_s.getfloat('Conductivity', 'mcond_c')}
        ndD_s["simPartCond"] = {"a": P_s.getboolean('Conductivity', 'simPartCond_a'),
                "c": P_s.getboolean('Conductivity', 'simPartCond_c')}
        dD_s["G_mean"] = {"a": P_s.getfloat('Conductivity', 'G_mean_a'),
                "c": P_s.getfloat('Conductivity', 'G_mean_c')}
        dD_s["G_stddev"] = {"a": P_s.getfloat('Conductivity', 'G_stddev_a'),
                "c": P_s.getfloat('Conductivity', 'G_stddev_c')}

        # Geometry
        dD_s["L"] = {"a": P_s.getfloat('Geometry', 'L_a'),
                "c": P_s.getfloat('Geometry', 'L_c'),
                "s": P_s.getfloat('Geometry', 'L_s')}
        ndD_s["P_L"] = {"a": P_s.getfloat('Geometry', 'P_L_a'),
                "c": P_s.getfloat('Geometry', 'P_L_c')}
        ndD_s["poros"] = {"a": P_s.getfloat('Geometry', 'poros_a'),
                "c": P_s.getfloat('Geometry', 'poros_c')}
        try:
            ndD_s["poros"]["s"] = P_s.getfloat('Geometry', 'poros_s')
        except ConfigParser.NoOptionError:
            ndD_s["poros"]["s"] = 1.

        # Electrolyte
        c0 = dD_s["c0"] = P_s.getfloat('Electrolyte', 'c0')
        zp = ndD_s["zp"] = P_s.getfloat('Electrolyte', 'zp')
        zm = ndD_s["zm"] = P_s.getfloat('Electrolyte', 'zm')
        if zm > 0:
            zm = ndD_s["zm"] = -zm
        nup = ndD_s["nup"] = P_s.getfloat('Electrolyte', 'nup')
        num = ndD_s["num"] = P_s.getfloat('Electrolyte', 'num')
        elyteModelType = ndD_s["elyteModelType"] = P_s.get('Electrolyte',
                'elyteModelType')
        SMset = ndD_s["SMset"] = P_s.get('Electrolyte', 'SMset')
        Dref = dD_s["Dref"] = elyte_CST.getProps(SMset)[-1]
        Dp = dD_s["Dp"] = P_s.getfloat('Electrolyte', 'Dp')
        Dm = dD_s["Dm"] = P_s.getfloat('Electrolyte', 'Dm')

        # Constants
        k = dD_s["k"] = 1.381e-23 # J/(K particle)
        Tref = dD_s["Tref"] = 298. # K
        e = dD_s["e"] = 1.602e-19 # C
        N_A = dD_s["N_A"] = 6.022e23 # particle/mol
        F = dD_s["F"] = dD_s["e"] * dD_s["N_A"] # C/mol

        for trode in ndD_s["trodes"]:
            dD = {}
            ndD = {}
            P = P_e[trode]

            Type = ndD["type"] = P.get('Particles', 'type')
            dD["disc"] = P.getfloat('Particles', 'discretization')
            Shape = ndD["shape"] = P.get('Particles', 'shape')
            dD["thickness"] = P.getfloat('Particles', 'thickness')
            if Type in ["ACR"]:
                ndD["simSurfCond"] = P.getboolean('Conductivity', 'simSurfCond')
                dD["scond"] = P.getfloat('Conductivity', 'scond')

            # Material
            # 1var and 2var parameters
#            if Type not in ["diffn"]:
            if Type in ["ACR", "CHR", "CHR2", "homog", "homog2",
                    "homog_sdn", "homog2_sdn"]:
                dD["Omga"] = P.getfloat('Material', 'Omega_a')
                dD["kappa"] = P.getfloat('Material', 'kappa')
                dD["B"] = P.getfloat('Material', 'B')
                dD["Vstd"] = P.getfloat('Material', 'Vstd')
            dD["rho_s"] = P.getfloat('Material', 'rho_s')
            ndD["delPhiEqFit"] = P.getboolean('Material', 'delPhiEqFit')
            if ndD["delPhiEqFit"]:
                ndD["delPhiFunc"] = P.get('Material', 'delPhiFunc')
            else:
                ndD["delPhiFunc"] = None
#            if Type not in ["ACR", "homog", "homog_sdn"]:
            if Type in ["diffn", "diffn2", "CHR", "CHR2"]:
                dD["Dsld"] = P.getfloat('Material', 'Dsld')
            if Type in ["CHR", "CHR2"]:
                dD["dgammadc"] = P.getfloat('Material', 'dgammadc')
            if Type in ["ACR"]:
                ndD["cwet"] = P.getfloat('Material', 'cwet')
#            if Type in ["diffn", "homog"]:
            ndD["logPad"] = P.getboolean('Material', 'logPad')
            ndD["noise"] = P.getboolean('Material', 'noise')

            # 2var (extra) parameters
            if Type in ["CHR2"]:
                dD["Omgb"] = P.getfloat('Material', 'Omega_b')
                dD["Omgc"] = P.getfloat('Material', 'Omega_c')
                dD["EvdW"] = P.getfloat('Material', 'EvdW')

            # Reactions
            ndD["rxnType"] = P.get('Reactions', 'rxnType')
            dD["k0"] = P.getfloat('Reactions', 'k0')
            ndD["alpha"] = P.getfloat('Reactions', 'alpha')
            dD["lambda"] = P.getfloat('Reactions', 'lambda')

            # electrode parameters
            dD_e[trode] = dD.copy()
            ndD_e[trode] = ndD.copy()

        # Post-processing
        self.test_system_input(dD_s, ndD_s)
        for trode in ndD_s["trodes"]:
            self.test_electrode_input(dD_e[trode], ndD_e[trode], dD_s, ndD_s)
        psd_raw, psd_num, psd_len, psd_area, psd_vol = self.distr_part(
                dD_s, ndD_s, dD_e, ndD_e)
        G = self.distr_G(dD_s, ndD_s)

        # Various calculated and defined parameters
        Lref = dD_s["Lref"] = dD_s["L"]["c"]
        # Ambipolar diffusivity
        Damb = dD_s["Damb"] = ((zp-zm)*Dp*Dm)/(zp*Dp-zm*Dm)
        # Cation transference number
        tp = ndD_s["tp"] = zp*Dp / (zp*Dp - zm*Dm)
        # Diffusive time scale
        if ndD_s["elyteModelType"] == "dilute":
            Dref = dD_s["Dref"] = Damb
        td = dD_s["td"] = Lref**2 / Dref

        # maximum concentration in electrode solids, mol/m^3
        # and electrode capacity ratio
        for trode in ndD_s['trodes']:
            dD_e[trode]['csmax'] = dD_e[trode]['rho_s'] / N_A # M
            dD_e[trode]['cap'] = (dD_s["e"] * dD_s['L'][trode] *
                    (1-ndD_s['poros'][trode]) *
                    ndD_s['P_L'][trode] *
                    dD_e[trode]['rho_s']) # C/m^2
        if "a" in ndD_s['trodes']:
            ndD_s['z'] = dD_e['c']['cap'] / dD_e['a']['cap']
        else:
            # flat plate anode with assumed infinite supply of metal
            ndD_s['z'] = 0.
        limtrode = ("c" if ndD_s["z"] < 1 else "a")
        CrateCurr = dD_s["CrateCurr"] = dD_e[limtrode]["cap"] / 3600. # A/m^2
        dD_s["currset"] = CrateCurr * dD_s["Crate"] # A/m^2

        # Some nondimensional parameters
        T = ndD_s["T"] = Tabs / Tref
        ndD_s["Rser"] = dD_s["Rser"] * e/(k*Tref) * (3600./td) * CrateCurr
        ndD_s["Dp"] = Dp / Dref
        ndD_s["Dm"] = Dm / Dref
        cref = dD_s["cref"] = 1000. # mol/m^3 = 1 M
        ndD_s["c0"] = c0 / cref
        ndD_s["phi_cathode"] = 0.
        ndD_s["currset"] = dD_s["Crate"]*td/3600
        ndD_s["Vset"] = dD_s["Vset"] * e/(k*Tref)
        ndD_s["tend"] = dD_s["tend"] / td
        if ndD_s["profileType"] == "CC" and not isClose(ndD_s["currset"], 0.):
            ndD_s["tend"] = ndD_s["capFrac"] / ndD_s["currset"]
        ndD_s["k0_foil"] = dD_s["k0_foil"] * (1./CrateCurr) * (td/3600.)

        # parameters which depend on the electrode
        dD_s["psd_raw"] = {}
        ndD_s["psd_num"] = {}
        dD_s["psd_len"] = {}
        dD_s["psd_area"] = {}
        dD_s["psd_vol"] = {}
        dD_s["G"] = {}
        ndD_s["G"] = {}
        ndD_s["psd_vol_FracVol"] = {}
        ndD_s["L"] = {}
        ndD_s["L"]["s"] = dD_s["L"]["s"] / Lref
        ndD_s["epsbeta"] = {}
        ndD_s["mcond"] = {}
        ndD_s["dphi_eq_ref"] = {}
#        ndD_s["phiRef"] = { # temporary, used for Vmax, Vmin
#                "a" : (e/(k*Tref))*dD_s["Vstd"]["a"],
#                "c" : 0.}
        ndD_s["phiRef"] = {"a" : 0.} # temporary, used for Vmax, Vmin
#        ndD_s["lambda"] = {}
#        ndD_s["B"] = {}
#        ndD_s["kappa"] = {}
#        ndD_s["k0"] = {}
#        ndD_s["beta_s"] = {}
#        ndD_s["delta_L"] = {}
#        ndD_s["scond"] = {}
#        ndD_s["Dsld"] = {}
#        ndD_s["Omga"] = {}
        for trode in ndD_s["trodes"]:
            ## System scale parameters
            # particle size distribution and connectivity distribution
            dD_s["psd_raw"][trode] = psd_raw[trode]
            ndD_s["psd_num"][trode] = psd_num[trode]
            dD_s["psd_len"][trode] = psd_len[trode]
            dD_s["psd_area"][trode] = psd_area[trode]
            dD_s["psd_vol"][trode] = psd_vol[trode]
            dD_s["G"][trode] = G[trode]
            # Sums of all particle volumes within each simulated
            # electrode volume
            Vuvec = np.sum(dD_s["psd_vol"][trode], axis=1)
            # Fraction of individual particle volume compared to total
            # volume of particles _within the simulated electrode
            # volume_
            ndD_s["psd_vol_FracVol"][trode] = (dD_s["psd_vol"][trode] /
                    Vuvec[:, np.newaxis])
            ndD_s["L"][trode] = dD_s["L"][trode]/Lref
            ndD_s["epsbeta"][trode] = (
                    (1-ndD_s['poros'][trode]) * ndD_s['P_L'][trode] *
                    dD_e[trode]["csmax"]/c0)
            ndD_s["mcond"][trode] = (
                    dD_s['mcond'][trode] * (td * k * N_A * Tref) /
                    (Lref**2 * F**2 * c0))
            vols = dD_s["psd_vol"][trode]
            ndD_s["G"][trode] = (dD_s["G"][trode] * (k*Tref/e) * td /
                    (F*dD_e[trode]["csmax"]*vols))

            ## Electrode particle parameters
            Type = ndD_e[trode]['type']
            if Type in ["diffn", "homog"] and ndD_e[trode]["delPhiEqFit"]:
                material = ndD_e[trode]['delPhiFunc']
                fits = delta_phi_fits.DPhiFits(ndD_s["T"])
                phifunc = fits.materialData[material]
                ndD_e[trode]["dphi_eq_ref"] = phifunc(ndD_s['cs0'][trode], 0)
                ndD_s["phiRef"][trode] = ndD_e[trode]["dphi_eq_ref"]
            else:
                ndD_e[trode]["dphi_eq_ref"] = 0.
                ndD_s["phiRef"][trode] = (e/(k*Tref))*dD_e[trode]["Vstd"]
            lmbda = ndD_e[trode]["lambda"] = dD_e[trode]["lambda"]/(k*Tref)
            if Type not in ["diffn"]:
                ndD_e[trode]["Omga"] = dD_e[trode]["Omga"] / (k*Tref)
                ndD_e[trode]["B"] = dD_e[trode]['B']/(k*Tref*dD_e[trode]['rho_s'])
            if Type in ["CHR2"]:
                ndD_e[trode]["Omgb"] = dD_e[trode]["Omgb"] / (k*Tref)
                ndD_e[trode]["Omgc"] = dD_e[trode]["Omgc"] / (k*Tref)
                ndD_e[trode]["EvdW"] = dD_e[trode]["EvdW"] / (k*Tref)
#            lens = dD["psd_len"][trode]
#            areas = dD["psd_area"][trode]
#            vols = dD["psd_vol"][trode]
            Nvol, Npart = psd_raw[trode].shape
            # Electrode parameters which depend on the individual
            # particle (e.g. because of non-dimensionalization)
            ndD_e[trode]["indvPart"] = np.empty((Nvol, Npart), dtype=object)
            for i in range(Nvol):
                for j in range(Npart):
                    # Bring in a copy of the nondimensional parameters
                    # we've calculated so far which are the same for
                    # all particles.
                    ndD_e[trode]["indvPart"][i, j] = ndD_e[trode].copy()
                    # This creates a reference for shorthand.
                    ndD_tmp = ndD_e[trode]["indvPart"][i, j]
                    # This specific particle dimensions
                    ndD_tmp["N"] = ndD_s["psd_num"][trode][i, j]
                    plen = dD_s["psd_len"][trode][i, j]
                    parea = dD_s["psd_area"][trode][i, j]
                    pvol = dD_s["psd_vol"][trode][i, j]
                    if Type not in ["diffn"]:
                        ndD_tmp["kappa"] = (dD_e[trode]['kappa'] /
                                (k*Tref*dD_e[trode]['rho_s']*plen**2))
                    if Type in ["CHR", "CHR2"]:
                        ndD_tmp["beta_s"] = (dD_e[trode]['dgammadc'] *
                                plen * dD_e[trode]['rho_s'] /
                                dD_e[trode]['kappa'])
                    if Type in ["ACR"]:
                        ndD_tmp["scond"] = (dD_e[trode]['scond'] * (k*Tref) /
                                (dD_e[trode]['k0']*e*plen**2))
                    if Type not in ["ACR", "homog", "homog_sdn"]:
                        ndD_tmp["Dsld"] = dD_e[trode]['Dsld']*td/plen**2
                    ndD_tmp["k0"] = ((parea/pvol)*dD_e[trode]['k0']*td
                            / (F*dD_e[trode]["csmax"]))
                    ndD_tmp["delta_L"] = pvol/(parea*plen)
                    if Type == "homog_sdn":
                        ndD_tmp["Omga"] = self.size2regsln(plen)
#            ndD["kappa"][trode] = (dD['kappa'][trode] /
#                    (k*Tref*dD['rho_s'][trode]*lens**2))
#            k0 = ndD["k0"][trode] = (
#                    ((areas/vols)*dD['k0'][trode]*td) /
#                    (F*dD["csmax"][trode]))
#            ndD["beta_s"][trode] = (dD['dgammasdc'][trode]*lens*
#                    dD['rho_s'][trode]/dD['kappa'][trode])
#            ndD["delta_L"][trode] = vols/(areas*lens)
#            ndD["scond"][trode] = (dD['scond'][trode] * (k*Tref) /
#                    (dD['k0'][trode]*e*lens**2))
#            ndD["Dsld"][trode] = dD['Dsld'][trode]*td/lens**2
#            solidType = ndD["solidType"][trode]
#            if solidType in ["homog", "ACR", "CHR", "diffn"]:
#                ndD["Omga"][trode] = (dD["Omga"][trode] / (k*Tref) *
#                        np.ones(psd_num[trode].shape))
#            elif solidType in ["homog_sdn"]:
#                # Not sure about factor of nondimensional T.
#                ndD["Omga"][trode] = T*self.size2regsln(lens)
#            else:
#                raise NotImplementedError("Solid types missing here")

        # Set up voltage cutoff values
        ndDVref = ndD_s["phiRef"]["c"] - ndD_s["phiRef"]["a"]
        ndD_s["phimin"] = -((e/(k*Tref))*dD_s["Vmax"] - ndDVref)
        ndD_s["phimax"] = -((e/(k*Tref))*dD_s["Vmin"] - ndDVref)

        return dD_s, ndD_s, dD_e, ndD_e

    def distr_part(self, dD_s, ndD_s, dD_e, ndD_e):
        psd_raw = {}
        psd_num = {}
        psd_len = {}
        psd_area = {}
        psd_vol = {}
        for trode in ndD_s["trodes"]:
            Nv = ndD_s["Nvol"][trode]
            Np = ndD_s["Npart"][trode]
            mean = dD_s["psd_mean"][trode]
            stddev = dD_s["psd_stddev"][trode]
            solidType = ndD_e[trode]["type"]
            # Make a length-sampled particle size distribution
            # Log-normally distributed
            if isClose(dD_s["psd_stddev"][trode], 0.):
                raw = (dD_s["psd_mean"][trode] *
                        np.ones((Nv, Np)))
            else:
                var = stddev**2
                mu = np.log((mean**2)/np.sqrt(var+mean**2))
                sigma = np.sqrt(np.log(var/(mean**2)+1))
                raw = np.random.lognormal(
                        mu, sigma, size=(Nv, Np))
            psd_raw[trode] = raw
            # For particles with internal profiles, convert psd to
            # integers -- number of steps
            solidDisc = dD_e[trode]["disc"]
            if solidType in ["ACR"]:
                psd_num[trode] = np.ceil(psd_raw[trode]/solidDisc).astype(np.integer)
                psd_len[trode] = solidDisc*psd_num[trode]
            elif solidType in ["CHR", "diffn", "CHR2", "diffn2"]:
                psd_num[trode] = np.ceil(psd_raw[trode]/solidDisc).astype(np.integer) + 1
                psd_len[trode] = solidDisc*(psd_num[trode] - 1)
            # For homogeneous particles (only one "volume" per particle)
            elif solidType in ["homog", "homog_sdn", "homog2", "homog2_sdn"]:
                # Each particle is only one volume
                psd_num[trode] = np.ones(psd_raw[trode].shape).astype(np.integer)
                # The lengths are given by the original length distr.
                psd_len[trode] = psd_raw[trode]
            else:
                raise NotImplementedError("Solid types missing here")
            # Calculate areas and volumes
            solidShape = ndD_e[trode]["shape"]
            if solidShape == "sphere":
                psd_area[trode] = (4*np.pi)*psd_len[trode]**2
                psd_vol[trode] = (4./3)*np.pi*psd_len[trode]**3
            elif solidShape == "C3":
                psd_area[trode] = 2 * 1.2263 * psd_len[trode]**2
                psd_vol[trode] = 1.2263 * psd_len[trode]**2 * dD_e[trode]['thickness']
            elif solidShape == "cylinder":
                psd_area[trode] = 2 * np.pi * psd_len[trode] * dD_e[trode]['thickness']
                psd_vol[trode] = np.pi * psd_len[trode]**2 * dD_e[trode]['thickness']
        return psd_raw, psd_num, psd_len, psd_area, psd_vol

    def distr_G(self, dD, ndD):
        G = {}
        for trode in ndD["trodes"]:
            Nv = ndD["Nvol"][trode]
            Np = ndD["Npart"][trode]
            mean = dD["G_mean"][trode]
            stddev = dD["G_stddev"][trode]
            if isClose(stddev, 0.):
                G[trode] = mean * np.ones((Nv, Np))
            else:
                var = stddev**2
                mu = np.log((mean**2)/np.sqrt(var+mean**2))
                sigma = np.sqrt(np.log(var/(mean**2)+1))
                G[trode] = np.random.lognormal(mu, sigma, size=(Nv, Np))
        return G

    def size2regsln(self, size):
        """
        This function returns the non-dimensional regular solution
        parameter which creates a barrier height that corresponds to
        the given particle size (C3 particle, measured in nm in the
        [100] direction). The barrier height vs size is taken from
        Cogswell and Bazant 2013, and the reg sln vs barrier height
        was done by TRF 2014 (Ferguson and Bazant 2014).
        """
        # First, this function wants the argument to be in [nm]
        size *= 1e+9
        # Parameters for polynomial curve fit
        p1 = -1.168e4
        p2 = 2985
        p3 = -208.3
        p4 = -8.491
        p5 = -10.25
        p6 = 4.516
        # The nucleation barrier depends on the ratio of the particle
        # wetted area to total particle volume.
        # *Wetted* area to volume ratio for C3 particles (Cogswell
        # 2013 or Kyle Smith)
        AV = 3.6338/size
        # Fit function (TRF, "SWCS" paper 2014)
        param = p1*AV**5 + p2*AV**4 + p3*AV**3 + p4*AV**2 + p5*AV + p6
        # replace values less than 2 with 2.
        if type(param) == np.ndarray:
            param[param < 2] = 2.
        else:
            if param < 2:
                param = 2.
        return param

    def test_system_input(self, dD, ndD):
        if not isClose(dD['Tabs'], 298.):
            raise Exception("Temp dependence not implemented")
        if ndD['Nvol']["c"] < 1:
            raise Exception("Must have at least one porous electrode")
        return

    def test_electrode_input(self, dD, ndD, dD_s, ndD_s):
        T298 = isClose(dD_s['Tabs'], 298.)
        solidType = ndD['type']
        solidShape = ndD['shape']
#        if ndD['simSurfCond'] and solidType != "ACR":
#            raise Exception("simSurfCond req. ACR")
        if solidType in ["ACR", "homog_sdn"] and solidShape != "C3":
            raise Exception("ACR and homog_sdn req. C3 shape")
        if (solidType in ["CHR", "diffn"] and solidShape not in
                ["sphere", "cylinder"]):
            raise NotImplementedError("CHR and diffn req. sphere or cylinder")
        if (solidType not in ndD_s["1varTypes"]) and (solidType not in
                ndD_s["2varTypes"]):
            raise NotImplementedError("Input solidType not defined")
        if solidShape not in ["C3", "sphere", "cylinder"]:
            raise NotImplementedError("Input solidShape not defined")
        if solidType == "homog_sdn" and not T298:
            raise NotImplementedError("homog_snd req. Tabs=298")
        if ndD["rxnType"] == "BV" and ndD_s["elyteModelType"] == "SM":
            raise Exception("BV currently requires dilute elyte model")
        try:
            if ndD['delPhiEqFit']:
                if ndD['delPhiFunc'] == "LiMn2O4" and not T298:
                    raise Exception("LiMn204 req. Tabs = 298 K")
                if ndD['delPhiFunc'] == "LiC6" and not T298:
                    raise Exception("LiC6 req. Tabs = 298 K")
                if ndD['delPhiFunc'] == "NCA1" and not T298:
                    raise Exception("NCA1 req. Tabs = 298 K")
                if solidType not in ["diffn", "homog"]:
                    raise NotImplementedError("delPhiEqFit req. solidType = diffn or homog")
        except KeyError:
            pass
        return

    def writeConfigFile(self, P, filename="input_params.cfg"):
        fo = open(filename, "w")
        P.write(fo)
        fo.close()
        return
    def writeDicts(self, dD, ndD, filenamebase="input_dict"):
        pickle.dump(dD, open(filenamebase + "_dD.p", "wb"))
        pickle.dump(ndD, open(filenamebase + "_ndD.p", "wb"))
        return
    def readDicts(self, filenamebase="input_dict"):
        dD = pickle.load(open(filenamebase + "_dD.p", "rb"))
        ndD = pickle.load(open(filenamebase + "_ndD.p", "rb"))
        return dD, ndD

def isClose(a, b):
    if np.abs(a - b) < 1e-12:
        return True
    else:
        return False
