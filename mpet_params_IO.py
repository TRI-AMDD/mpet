import ConfigParser
import pickle

import numpy as np
import scipy.special as spcl

import delta_phi_fits

class mpetIO():

    def getConfig(self, paramfile="params.cfg"):
        P = ConfigParser.RawConfigParser()
        P.optionxform = str
        P.read(paramfile)
        return P

    def getDictFromConfig(self, P):
        # Dictionary of dimensional parameters
        dD = {}
        # Dictionary of nondimensional parameters
        ndD = {}

        # Simulation Parameters
        ndD["profileType"] = P.get('Sim Params', 'profileType')
        dD["Crate"] = P.getfloat('Sim Params', 'Crate')
        dD["Vset"] = P.getfloat('Sim Params', 'Vset')
        ndD["capFrac"] = P.getfloat('Sim Params', 'capFrac')
        dD["tend"] = P.getfloat('Sim Params', 'tend')
        ndD["tsteps"] = P.getfloat('Sim Params', 'tsteps')
        Tabs = dD["Tabs"] = P.getfloat('Sim Params', 'T')
        ndD["Nvol"] = {"a": P.getint('Sim Params', 'Nvol_a'),
                       "c": P.getint('Sim Params', 'Nvol_c'),
                       "s": P.getint('Sim Params', 'Nvol_s')}
        ndD["Npart"] = {"a": P.getint('Sim Params', 'Npart_a'),
                        "c": P.getint('Sim Params', 'Npart_c')}
        ndD["trodes"] = ["c"]
        if ndD["Nvol"]["a"] >= 1:
            ndD["trodes"].append("a")

        # Particle info
        dD["psd_mean"] = {"a": P.getfloat('Particles', 'mean_a'),
                          "c":  P.getfloat('Particles', 'mean_c')}
        dD["psd_stddev"] = {"a": P.getfloat('Particles', 'stddev_a'),
                            "c": P.getfloat('Particles', 'stddev_c')}
        ndD["cs0"] = {"a": P.getfloat('Particles', 'cs0_a'),
                "c": P.getfloat('Particles', 'cs0_c')}
        ndD["solidType"] = {"a": P.get('Particles', 'solidType_a'),
                "c": P.get('Particles', 'solidType_c')}
        dD["solidDisc"] = {"a": P.getfloat('Particles', 'solidDisc_a'),
                "c": P.getfloat('Particles', 'solidDisc_c')}
        ndD["solidShape"] = {"a": P.get('Particles', 'solidShape_a'),
                "c": P.get('Particles', 'solidShape_c')}
        dD["partThick"] = {"a": P.getfloat('Particles', 'partThick_a'),
                "c": P.getfloat('Particles', 'partThick_c')}

        # Conductivity
        ndD["simBulkCond"] = {"a": P.getboolean('Conductivity', 'simBulkCond_a'),
                "c": P.getboolean('Conductivity', 'simBulkCond_c')}
        dD["mcond"] = {"a": P.getfloat('Conductivity', 'mcond_a'),
                "c": P.getfloat('Conductivity', 'mcond_c')}
        ndD["simPartCond"] = {"a": P.getboolean('Conductivity', 'simPartCond_a'),
                "c": P.getboolean('Conductivity', 'simPartCond_c')}
        dD["G_mean"] = {"a": P.getfloat('Conductivity', 'G_mean_a'),
                "c": P.getfloat('Conductivity', 'G_mean_c')}
        dD["G_stddev"] = {"a": P.getfloat('Conductivity', 'G_stddev_a'),
                "c": P.getfloat('Conductivity', 'G_stddev_c')}
        ndD["simSurfCond"] = {"a": P.getboolean('Conductivity', 'simSurfCond_a'),
                "c": P.getboolean('Conductivity', 'simSurfCond_c')}
        dD["scond"] = {"a": P.getfloat('Conductivity', 'scond_a'),
                "c": P.getfloat('Conductivity', 'scond_c')}

        # Materials
        dD["Omga"] = {"a": P.getfloat('Materials', 'Omega_a_a'),
                "c": P.getfloat('Materials', 'Omega_a_c')}
        dD["kappa"] = {"a": P.getfloat('Materials', 'kappa_a'),
                "c": P.getfloat('Materials', 'kappa_c')}
        dD["B"] = {"a": P.getfloat('Materials', 'B_a'),
                "c": P.getfloat('Materials', 'B_c')}
        dD["rhos"] = {"a": P.getfloat('Materials', 'rhos_a'),
                "c": P.getfloat('Materials', 'rhos_c')}
        dD["Vstd"] = {"a": P.getfloat('Materials', 'Vstd_a'),
                "c": P.getfloat('Materials', 'Vstd_c')}
        dD["Dsld"] = {"a": P.getfloat('Materials', 'Dsld_a'),
                "c": P.getfloat('Materials', 'Dsld_c')}
        dD["dgammasdc"] = {"a": P.getfloat('Materials', 'dgammasdc_a'),
                "c": P.getfloat('Materials', 'dgammasdc_c')}
        ndD["cwet"] = {"a": P.getfloat('Materials', 'cwet_a'),
                "c": P.getfloat('Materials', 'cwet_c')}
        ndD["delPhiEqFit"] = {"a": P.getboolean('Materials', 'delPhiEqFit_a'),
                "c": P.getboolean('Materials', 'delPhiEqFit_c')}
        ndD["material"] = {"a": P.get('Materials', 'material_a'),
                "c": P.get('Materials', 'material_c')}

        # Reactions
        ndD["rxnType"] = {"a": P.get('Reactions', 'rxnType_a'),
                "c": P.get('Reactions', 'rxnType_c')}
        dD["k0"] = {"a": P.getfloat('Reactions', 'k0_a'),
                "c": P.getfloat('Reactions', 'k0_c')}
        ndD["alpha"] = {"a": P.getfloat('Reactions', 'alpha_a'),
                "c": P.getfloat('Reactions', 'alpha_c')}
        dD["lambda"] = {"a": P.getfloat('Reactions', 'lambda_a'),
                "c": P.getfloat('Reactions', 'lambda_c')}

        # Geometry
        dD["L"] = {"a": P.getfloat('Geometry', 'L_a'),
                "c": P.getfloat('Geometry', 'L_c'),
                "s": P.getfloat('Geometry', 'L_s')}
        ndD["P_L"] = {"a": P.getfloat('Geometry', 'P_L_a'),
                "c": P.getfloat('Geometry', 'P_L_c')}
        ndD["poros"] = {"a": P.getfloat('Geometry', 'poros_a'),
                "c": P.getfloat('Geometry', 'poros_c')}
        ndD["poros"]["s"] = 1.

        # Electrolyte
        c0 = dD["c0"] = P.getfloat('Electrolyte', 'c0')
        Dp = dD["Dp"] = P.getfloat('Electrolyte', 'Dp')
        Dm = dD["Dm"] = P.getfloat('Electrolyte', 'Dm')
        zp = ndD["zp"] = P.getfloat('Electrolyte', 'zp')
        zm = ndD["zm"] = P.getfloat('Electrolyte', 'zm')

        # Constants
        k = dD["k"] = P.getfloat('Constants', 'k')
        Tref = dD["Tref"] = P.getfloat('Constants', 'Tref')
        e = dD["e"] = P.getfloat('Constants', 'e')
        N_A = dD["N_A"] = P.getfloat('Constants', 'N_A')
        F = dD["F"] = dD["e"] * dD["N_A"]

        # Post-processing
        self.test_input(dD, ndD)
        psd_raw, psd_num, psd_len, psd_area, psd_vol = self.distr_part(dD, ndD)
        G = self.distr_G(dD, ndD)

        # Various calculated and defined parameters
        Lref = dD["Lref"] = dD["L"]["c"]
        # maximum concentration in electrode solids, mol/m^3
        dD["csmax"] = {"a": dD['rhos']["a"]/N_A,
                "c": dD['rhos']["c"]/N_A}
        # Ambipolar diffusivity
        Damb = dD["Damb"] = ((zp+zm)*Dp*Dm)/(zp*Dp+zm*Dm)
        # Cation transference number
        tp = ndD["tp"] = zp*Dp / (zp*Dp + zm*Dm)
        # Diffusive time scale
        td = dD["td"] = Lref**2 / Damb
        # Electrode capacity ratio
        dD["cap"] = {}
        dD["cap"]["c"] = (dD['L']["c"] * (1-ndD['poros']["c"]) *
                ndD['P_L']["c"] * dD['rhos']["c"])
        if "a" in ndD["trodes"]:
            # full porous anode with finite capacity
            dD["cap"]["a"] = (dD['L']["a"] * (1-ndD['poros']["a"]) *
                    ndD['P_L']["a"] * dD['rhos']["a"])
            ndD["z"] = dD["cap"]["c"] / dD["cap"]["a"]
        else:
            # flat plate anode with assumed infinite supply of metal
            ndD["z"] = 0

        # Some nondimensional parameters
        T = ndD["T"] = Tabs/Tref
        ndD["Dp"] = Dp / Damb
        ndD["Dm"] = Dm / Damb
        ndD["c0"] = c0 / 1000. # normalize by 1 M
        ndD["phi_cathode"] = 0.
        ndD["currset"] = dD["Crate"]*td/3600
        ndD["Vset"] = dD["Vset"] * e/(k*Tref)
        ndD["tend"] = dD["tend"] / td
        if ndD["profileType"] == "CC" and ndD["currset"] != 0.0:
            ndD["tend"] = ndD["capFrac"] / ndD["currset"]
        else: # CV or zero current simulation
            ndD["tend"] = dD["tend"] / td

        # parameters which depend on the electrode
        dD["psd_raw"] = {}
        ndD["psd_num"] = {}
        dD["psd_len"] = {}
        dD["psd_area"] = {}
        dD["psd_vol"] = {}
        dD["G"] = {}
        ndD["psd_vol_FracTot"] = {}
        ndD["psd_vol_FracVol"] = {}
        ndD["L"] = {}
        ndD["L"]["s"] = dD["L"]["s"] / Lref
        ndD["epsbeta"] = {}
        ndD["mcond"] = {}
        ndD["dphi_eq_ref"] = {}
        ndD["lambda"] = {}
        ndD["MHC_erfstretch"] = {}
        ndD["B"] = {}
        ndD["kappa"] = {}
        ndD["k0"] = {}
        ndD["beta_s"] = {}
        ndD["delta_L"] = {}
        ndD["MHC_Aa"] = {}
        ndD["scond"] = {}
        ndD["Dsld"] = {}
        ndD["G"] = {}
        ndD["Omga"] = {}
        for trode in ndD["trodes"]:
            # particle size distribution and connectivity distribution
            dD["psd_raw"][trode] = psd_raw[trode]
            ndD["psd_num"][trode] = psd_num[trode]
            dD["psd_len"][trode] = psd_len[trode]
            dD["psd_area"][trode] = psd_area[trode]
            dD["psd_vol"][trode] = psd_vol[trode]
            dD["G"][trode] = G[trode]
            # Fraction of individual particle volume compared to total
            # particle volume
            ndD["psd_vol_FracTot"][trode] = (dD["psd_vol"][trode] /
                    np.sum(dD["psd_vol"][trode]))
            # Sums of all particle volumes within each simulated
            # electrode volume
            Vuvec = np.sum(dD["psd_vol"][trode], axis=1)
            # Fraction of individual particle volume compared to total
            # volume of particles _within the simulated electrode
            # volume_
            ndD["psd_vol_FracVol"][trode] = (dD["psd_vol"][trode] /
                    Vuvec[:, np.newaxis])
            ndD["L"][trode] = dD["L"][trode]/Lref
            ndD["epsbeta"][trode] = (
                    (1-ndD['poros'][trode]) * ndD['P_L'][trode] *
                    dD["csmax"][trode]/c0)
            ndD["mcond"][trode] = (
                    dD['mcond'][trode] * (td * k * N_A * Tref) /
                    (Lref**2 * F**2 * c0))
            if ndD["delPhiEqFit"][trode]:
                material = ndD['material'][trode]
                fits = delta_phi_fits.DPhiFits(ndD["T"])
                phifunc = fits.materialData[material]
                ndD["dphi_eq_ref"][trode] = phifunc(dD['cs0'][trode], 0)
            else:
                ndD["dphi_eq_ref"][trode] = 0.
            lmbda = ndD["lambda"][trode] = dD["lambda"][trode]/(k*Tref)
            MHC_erf_b = ndD["MHC_erfstretch"][trode] = 2*np.sqrt(lmbda)
            ndD["B"][trode] = dD['B'][trode]/(k*Tref*dD['rhos'][trode])
            lens = dD["psd_len"][trode]
            areas = dD["psd_area"][trode]
            vols = dD["psd_vol"][trode]
            ndD["kappa"][trode] = (dD['kappa'][trode] /
                    (k*Tref*dD['rhos'][trode]*lens**2))
            k0 = ndD["k0"][trode] = (
                    ((areas/vols)*dD['k0'][trode]*td) /
                    (F*dD["csmax"][trode]))
            ndD["beta_s"][trode] = (dD['dgammasdc'][trode]*lens*
                    dD['rhos'][trode]/dD['kappa'][trode])
            ndD["delta_L"][trode] = vols/(areas*lens)
            ndD["MHC_Aa"][trode] = k0 / (spcl.erf(-lmbda/MHC_erf_b) + 1)
            ndD["scond"][trode] = (dD['scond'][trode] * (k*Tref) /
                    (dD['k0'][trode]*e*lens**2))
            ndD["Dsld"][trode] = dD['Dsld'][trode]*td/lens**2
            ndD["G"][trode] = (dD["G"][trode] * (k*Tref/e) * td /
                    (F*dD["csmax"][trode]*vols))
            solidType = ndD["solidType"][trode]
            if solidType in ["homog", "ACR", "CHR", "diffn"]:
                ndD["Omga"][trode] = (dD["Omga"][trode] / (k*Tref) *
                        np.ones(psd_num[trode].shape))
            elif solidType in ["homog_sdn"]:
                # Not sure about factor of nondimensional T.
                ndD["Omga"][trode] = T*self.size2regsln(lens)
            else:
                raise NotImplementedError("Solid types missing here")

        return dD, ndD

    def distr_part(self, dD, ndD):
        psd_raw = {}
        psd_num = {}
        psd_len = {}
        psd_area = {}
        psd_vol = {}
        for trode in ndD["trodes"]:
            Nv = ndD["Nvol"][trode]
            Np = ndD["Npart"][trode]
            mean = dD["psd_mean"][trode]
            stddev = dD["psd_stddev"][trode]
            solidType = ndD["solidType"][trode]
            # Make a length-sampled particle size distribution
            # Log-normally distributed
            if dD["psd_mean"][trode] == 0:
                raw = (dD["psd_mean"][trode] *
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
            solidDisc = dD["solidDisc"][trode]
            if solidType in ["ACR"]:
                psd_num[trode] = np.ceil(psd_raw[trode]/solidDisc).astype(np.integer)
                psd_len[trode] = solidDisc*psd_num[trode]
            elif solidType in ["CHR", "diffn"]:
                psd_num[trode] = np.ceil(psd_raw[trode]/solidDisc).astype(np.integer) + 1
                psd_len[trode] = solidDisc*(psd_num[trode] - 1)
            # For homogeneous particles (only one "volume" per particle)
            elif solidType in ["homog", "homog_sdn"]:
                # Each particle is only one volume
                psd_num[trode] = np.ones(psd_raw[trode].shape).astype(np.integer)
                # The lengths are given by the original length distr.
                psd_len[trode] = psd_raw[trode]
            else:
                raise NotImplementedError("Solid types missing here")
            # Calculate areas and volumes
            solidShape = ndD["solidShape"][trode]
            if solidShape == "sphere":
                psd_area[trode] = (4*np.pi)*psd_len[trode]**2
                psd_vol[trode] = (4./3)*np.pi*psd_len[trode]**3
            elif solidShape == "C3":
                psd_area[trode] = 2 * 1.2263 * psd_len[trode]**2
                psd_vol[trode] = 1.2263 * psd_len[trode]**2 * dD['partThick'][trode]
            elif solidShape == "cylinder":
                psd_area[trode] = 2 * np.pi * psd_len[trode] * dD['partThick'][trode]
                psd_vol[trode] = np.pi * psd_len[trode]**2 * dD['partThick'][trode]
        return psd_raw, psd_num, psd_len, psd_area, psd_vol

    def distr_G(self, dD, ndD):
        G = {}
        for trode in ndD["trodes"]:
            Nv = ndD["Nvol"][trode]
            Np = ndD["Npart"][trode]
            mean = dD["G_mean"][trode]
            stddev = dD["G_stddev"][trode]
            if stddev == 0:
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
        Cogswell 2013, and the reg sln vs barrier height was done by
        TRF 2014.
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
        param[param < 2] = 2
        return param

    def test_input(self, dD, ndD):
        if dD['Tabs'] != 298 or dD['Tref'] != 298:
            raise Exception("Temp dependence not implemented")
        if ndD['Nvol']["c"] < 1:
            raise Exception("Must have at least one porous electrode")
        for trode in ndD["trodes"]:
            solidType = ndD['solidType'][trode]
            solidShape = ndD['solidShape'][trode]
            if ndD['simSurfCond'][trode] and solidType != "ACR":
                raise Exception("simSurfCond req. ACR")
            if solidType in ["ACR", "homog_sdn"] and solidShape != "C3":
                raise Exception("ACR and homog_sdn req. C3 shape")
            if solidType in ["CHR"] and solidShape not in ["sphere", "cylinder"]:
                raise NotImplementedError("CHR req. sphere or cylinder")
            if solidType not in ["ACR", "CHR", "homog", "homog_sdn", "diffn"]:
                raise NotImplementedError("Input solidType not defined")
            if solidShape not in ["C3", "sphere", "cylinder"]:
                raise NotImplementedError("Input solidShape not defined")
            if solidType == "homog_sdn" and (dD['Tabs'] != 298 or
                    dD['Tref'] != 298):
                raise NotImplementedError("homog_snd req. Tref=Tabs=298")
            if solidType in ["diffn"] and solidShape != "sphere":
                raise NotImplementedError("diffn currently req. sphere")
            if ndD['delPhiEqFit'][trode] and solidType not in ["diffn", "homog"]:
                if ndD['material'][trode] == "LiMn2O4" and dD['Tabs'] != 298:
                    raise Exception("LiMn204 req. Tabs = 298 K")
                if ndD['material'][trode] == "LiC6" and dD['Tabs'] != 298:
                    raise Exception("LiC6 req. Tabs = 298 K")
                if ndD['material'][trode] == "NCA1" and dD['Tabs'] != 298:
                    raise Exception("NCA1 req. Tabs = 298 K")
                raise NotImplementedError("delPhiEqFit req. solidType = diffn or homog")
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
