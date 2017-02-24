"""Helper functions for interacting with system parameters.

This module provides functions for various data format exchanges:
 - config files (on disk) <--> dictionaries of parameters (in memory)
 - dictionaries of parameters (in memory) <--> dictionaries of parameters ("pickled" on disk)

It also has various other functions used in processes for things such as generating
distributions from input means and standard deviations.
"""
import ast
import configparser
import os
import pickle

import numpy as np

import mpet.props_am as props_am
import mpet.props_elyte as props_elyte


def getConfigs(paramfile="params.cfg"):
    # system-level config
    if not os.path.isfile(paramfile):
        raise Exception("System param file doesn't exist!")
    P_s = getConfig(paramfile)

    # electrode config(s)
    P_e = {}
    # cathode config
    paramfile_c = P_s.get('Electrodes', 'cathode')
    # Look in the same directory as the input paramfile, wherever
    # that is.
    paramfileLoc = os.path.split(paramfile)[0]
    paramfile_c = os.path.join(paramfileLoc, paramfile_c)
    if not os.path.isfile(paramfile_c):
        raise Exception("Cathode param file doesn't exist!")
    P_e["c"] = getConfig(paramfile_c)

    # anode config
    if P_s.getint('Sim Params', 'Nvol_a') >= 1:
        paramfile_a = P_s.get('Electrodes', 'anode')
        paramfile_a = os.path.join(paramfileLoc, paramfile_a)
        if not os.path.isfile(paramfile_a):
            raise Exception("Anode param file doesn't exist!")
        P_e["a"] = getConfig(paramfile_a)
    return P_s, P_e


def getDictsFromConfigs(P_s, P_e):
    # Dictionary of dimensional system and electrode parameters
    dD_s = {}
    dD_e = {}
    # Dictionary of nondimensional system and electrode parameters
    ndD_s = {}
    ndD_e = {}

    # General particle classification
    ndD_s["2varTypes"] = ["diffn2", "CHR2", "homog2", "homog2_sdn"]
    ndD_s["1varTypes"] = ["ACR", "diffn", "CHR", "homog", "homog_sdn"]

    # Simulation parameters
    ndD_s["profileType"] = P_s.get('Sim Params', 'profileType')
    dD_s["Crate"] = P_s.getfloat('Sim Params', 'Crate')
    segs = dD_s["segments"] = ast.literal_eval(
        P_s.get('Sim Params', 'segments'))
    ndD_s["tramp"] = dD_s["tramp"] = P_s.getfloat('Sim Params', 'tramp')
    numsegs = dD_s["numsegments"] = len(segs)
    dD_s["Vmax"] = P_s.getfloat('Sim Params', 'Vmax')
    dD_s["Vmin"] = P_s.getfloat('Sim Params', 'Vmin')
    # Should have deprecation warnings to encourage users to
    # update their params files to mirror options in
    # configDefaults.
    ndD_s["relTol"] = P_s.getfloat('Sim Params', 'relTol')
    ndD_s["absTol"] = P_s.getfloat('Sim Params', 'absTol')
    ndD_s["randomSeed"] = P_s.getboolean('Sim Params', 'randomSeed')
    ndD_s["seed"] = P_s.getint('Sim Params', 'seed')
    if ndD_s["randomSeed"]:
        # This should affect all calls to np.random throughout the
        # simulation.
        np.random.seed(ndD_s["seed"])
    dD_s["Vset"] = P_s.getfloat('Sim Params', 'Vset')
    ndD_s["capFrac"] = P_s.getfloat('Sim Params', 'capFrac')
    dD_s["tend"] = P_s.getfloat('Sim Params', 'tend')
    ndD_s["prevDir"] = P_s.get('Sim Params', 'prevDir')
    ndD_s["tsteps"] = P_s.getfloat('Sim Params', 'tsteps')
    Tabs = dD_s["Tabs"] = P_s.getfloat('Sim Params', 'T')
    dD_s["Rser"] = P_s.getfloat('Sim Params', 'Rser')
    ndD_s["Nvol"] = {"a": P_s.getint('Sim Params', 'Nvol_a'),
                     "c": P_s.getint('Sim Params', 'Nvol_c'),
                     "s": P_s.getint('Sim Params', 'Nvol_s')}
    ndD_s["Npart"] = {"a": P_s.getint('Sim Params', 'Npart_a'),
                      "c": P_s.getint('Sim Params', 'Npart_c')}
    ndD_s["trodes"] = ["c"]
    if ndD_s["Nvol"]["a"] >= 1:
        ndD_s["trodes"].append("a")
    dD_s["k0_foil"] = P_s.getfloat('Electrodes', 'k0_foil')
    dD_s["Rfilm_foil"] = P_s.getfloat('Electrodes', 'Rfilm_foil')

    # Particle info
    dD_s["psd_mean"] = {"a": P_s.getfloat('Particles', 'mean_a'),
                        "c": P_s.getfloat('Particles', 'mean_c')}
    dD_s["psd_stddev"] = {"a": P_s.getfloat('Particles', 'stddev_a'),
                          "c": P_s.getfloat('Particles', 'stddev_c')}
    ndD_s["cs0"] = {"a": P_s.getfloat('Particles', 'cs0_a'),
                    "c": P_s.getfloat('Particles', 'cs0_c')}

    # Conductivity
    ndD_s["simBulkCond"] = {
        "a": P_s.getboolean('Conductivity', 'simBulkCond_a'),
        "c": P_s.getboolean('Conductivity', 'simBulkCond_c')}
    dD_s["mcond"] = {"a": P_s.getfloat('Conductivity', 'mcond_a'),
                     "c": P_s.getfloat('Conductivity', 'mcond_c')}
    ndD_s["simPartCond"] = {
        "a": P_s.getboolean('Conductivity', 'simPartCond_a'),
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
    ndD_s["poros"]["s"] = P_s.getfloat('Geometry', 'poros_s')
    ndD_s["BruggExp"] = {"a": P_s.getfloat('Geometry', 'BruggExp_a'),
                         "c": P_s.getfloat('Geometry', 'BruggExp_c'),
                         "s": P_s.getfloat('Geometry', 'BruggExp_s')}

    # Electrolyte
    c0 = dD_s["c0"] = P_s.getfloat('Electrolyte', 'c0')
    zp = ndD_s["zp"] = P_s.getfloat('Electrolyte', 'zp')
    zm = ndD_s["zm"] = P_s.getfloat('Electrolyte', 'zm')
    if zm > 0:
        zm = ndD_s["zm"] = -zm
    ndD_s["nup"] = P_s.getfloat('Electrolyte', 'nup')
    ndD_s["num"] = P_s.getfloat('Electrolyte', 'num')
    ndD_s["elyteModelType"] = P_s.get('Electrolyte', 'elyteModelType')
    SMset = ndD_s["SMset"] = P_s.get('Electrolyte', 'SMset')
    D_ref = dD_s["D_ref"] = dD_s["Dref"] = props_elyte.getProps(SMset)[-1]
    ndD_s["n_refTrode"] = P_s.getfloat('Electrolyte', 'n')
    ndD_s["sp"] = P_s.getfloat('Electrolyte', 'sp')
    Dp = dD_s["Dp"] = P_s.getfloat('Electrolyte', 'Dp')
    Dm = dD_s["Dm"] = P_s.getfloat('Electrolyte', 'Dm')

    # Constants
    k = dD_s["k"] = 1.381e-23  # J/(K particle)
    T_ref = dD_s["T_ref"] = dD_s["Tref"] = 298.  # K
    e = dD_s["e"] = 1.602e-19  # C
    N_A = dD_s["N_A"] = 6.022e23  # particle/mol
    F = dD_s["F"] = dD_s["e"] * dD_s["N_A"]  # C/mol

    for trode in ndD_s["trodes"]:
        dD = {}
        ndD = {}
        P = P_e[trode]

        # Particles
        Type = ndD["type"] = P.get('Particles', 'type')
        dD["disc"] = P.getfloat('Particles', 'discretization')
        ndD["shape"] = P.get('Particles', 'shape')
        dD["thickness"] = P.getfloat('Particles', 'thickness')

        # Material
        # both 1var and 2var parameters
        dD["Omga"] = P.getfloat('Material', 'Omega_a')
        dD["Omgb"] = P.getfloat('Material', 'Omega_b')
        dD["Omgc"] = P.getfloat('Material', 'Omega_c')
        dD["kappa"] = P.getfloat('Material', 'kappa')
        dD["B"] = P.getfloat('Material', 'B')
        dD["EvdW"] = P.getfloat('Material', 'EvdW')
        dD["rho_s"] = P.getfloat('Material', 'rho_s')
        dD["D"] = P.getfloat('Material', 'D')
        ndD["Dfunc"] = P.get('Material', 'Dfunc')
        dD["dgammadc"] = P.getfloat('Material', 'dgammadc')
        ndD["cwet"] = P.getfloat('Material', 'cwet')
        ndD["muRfunc"] = P.get('Material', 'muRfunc')
        ndD["logPad"] = P.getboolean('Material', 'logPad')
        ndD["noise"] = P.getboolean('Material', 'noise')
        ndD["noise_prefac"] = P.getfloat('Material', 'noise_prefac')
        ndD["numnoise"] = P.getint('Material', 'numnoise')

        # Reactions
        ndD["rxnType"] = P.get('Reactions', 'rxnType')
        dD["k0"] = P.getfloat('Reactions', 'k0')
        ndD["alpha"] = P.getfloat('Reactions', 'alpha')
        dD["lambda"] = P.getfloat('Reactions', 'lambda')
        dD["Rfilm"] = P.getfloat('Reactions', 'Rfilm')

        # electrode parameters
        dD_e[trode] = dD.copy()
        ndD_e[trode] = ndD.copy()

    # Post-processing
    test_system_input(dD_s, ndD_s)
    for trode in ndD_s["trodes"]:
        test_electrode_input(dD_e[trode], ndD_e[trode], dD_s, ndD_s)
    psd_raw, psd_num, psd_len, psd_area, psd_vol = distr_part(
        dD_s, ndD_s, dD_e, ndD_e)
    G = distr_G(dD_s, ndD_s)

    # Various calculated and defined parameters
    L_ref = dD_s["L_ref"] = dD_s["Lref"] = dD_s["L"]["c"]
    c_ref = dD_s["c_ref"] = dD_s["cref"] = 1000.  # mol/m^3 = 1 M
    # Ambipolar diffusivity
    Damb = dD_s["Damb"] = ((zp-zm)*Dp*Dm)/(zp*Dp-zm*Dm)
    # Cation transference number
    ndD_s["tp"] = zp*Dp / (zp*Dp - zm*Dm)
    # Diffusive time scale
    if ndD_s["elyteModelType"] == "dilute":
        D_ref = dD_s["D_ref"] = Damb
    t_ref = dD_s["t_ref"] = dD_s["td"] = L_ref**2 / D_ref
    curr_ref = dD_s["curr_ref"] = 3600. / t_ref
    dD_s["mcond_ref"] = (L_ref**2 * F**2 * c_ref) / (t_ref * k * N_A * T_ref)
    dD_s["elytei_ref"] = F*c_ref*D_ref / L_ref
    # maximum concentration in electrode solids, mol/m^3
    # and electrode capacity ratio
    for trode in ndD_s['trodes']:
        dD_e[trode]['csmax'] = dD_e[trode]['rho_s'] / N_A
        if ndD_e[trode]['type'] in ndD_s['1varTypes']:
            dD_e[trode]['cs_ref'] = dD_e[trode]['csmax']
        elif ndD_e[trode]['type'] in ndD_s['2varTypes']:
            dD_e[trode]['cs_ref'] = 0.5 * dD_e[trode]['csmax']
        dD_e[trode]['cap'] = (
            dD_s["e"] * dD_s['L'][trode] * (1-ndD_s['poros'][trode])
            * ndD_s['P_L'][trode] * dD_e[trode]['rho_s'])  # C/m^2
    if "a" in ndD_s['trodes']:
        ndD_s['z'] = dD_e['c']['cap'] / dD_e['a']['cap']
    else:
        # flat plate anode with assumed infinite supply of metal
        ndD_s['z'] = 0.
    limtrode = ("c" if ndD_s["z"] < 1 else "a")
    CrateCurr = dD_s["CrateCurr"] = dD_e[limtrode]["cap"] / 3600.  # A/m^2
    dD_s["currset"] = CrateCurr * dD_s["Crate"]  # A/m^2
    Rser_ref = dD_s["Rser_ref"] = (k*T_ref/e) / (curr_ref*CrateCurr)

    # Some nondimensional parameters
    ndD_s["T"] = Tabs / T_ref
    ndD_s["Rser"] = dD_s["Rser"] / Rser_ref
    ndD_s["Dp"] = Dp / D_ref
    ndD_s["Dm"] = Dm / D_ref
    ndD_s["c0"] = c0 / c_ref
    ndD_s["phi_cathode"] = 0.
    ndD_s["currset"] = dD_s["Crate"] / curr_ref
    ndD_s["k0_foil"] = dD_s["k0_foil"] * (1./(curr_ref*CrateCurr))
    ndD_s["Rfilm_foil"] = dD_s["Rfilm_foil"] / Rser_ref

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
    ndD_s["L"]["s"] = dD_s["L"]["s"] / L_ref
    ndD_s["beta"] = {}
    ndD_s["mcond"] = {}
    ndD_s["phiRef"] = {"a": 0.}  # temporary, used for Vmax, Vmin
    for trode in ndD_s["trodes"]:
        # System scale parameters
        # particle size distribution and connectivity distribution
        if ndD_s["prevDir"] == "false":
            dD_s["psd_raw"][trode] = psd_raw[trode]
            ndD_s["psd_num"][trode] = psd_num[trode]
            dD_s["psd_len"][trode] = psd_len[trode]
            dD_s["psd_area"][trode] = psd_area[trode]
            dD_s["psd_vol"][trode] = psd_vol[trode]
            dD_s["G"][trode] = G[trode]
        else:
            dD_sPrev, ndD_sPrev = readDicts(
                os.path.join(ndD_s["prevDir"], "input_dict_system"))
            dD_s["psd_raw"][trode] = dD_sPrev["psd_raw"][trode]
            ndD_s["psd_num"][trode] = ndD_sPrev["psd_num"][trode]
            dD_s["psd_len"][trode] = dD_sPrev["psd_len"][trode]
            dD_s["psd_area"][trode] = dD_sPrev["psd_area"][trode]
            dD_s["psd_vol"][trode] = dD_sPrev["psd_vol"][trode]
            dD_s["G"][trode] = dD_sPrev["G"][trode]
        # Sums of all particle volumes within each simulated
        # electrode volume
        Vuvec = np.sum(dD_s["psd_vol"][trode], axis=1)
        # Fraction of individual particle volume compared to total
        # volume of particles _within the simulated electrode
        # volume_
        ndD_s["psd_vol_FracVol"][trode] = (
            dD_s["psd_vol"][trode] / Vuvec[:, np.newaxis])
        ndD_s["L"][trode] = dD_s["L"][trode]/L_ref
        ndD_s["beta"][trode] = dD_e[trode]["csmax"]/c_ref
        ndD_s["mcond"][trode] = dD_s['mcond'][trode] / dD_s["mcond_ref"]
        vols = dD_s["psd_vol"][trode]
        ndD_s["G"][trode] = (
            dD_s["G"][trode] * (k*T_ref/e) * t_ref
            / (F*dD_e[trode]["csmax"]*vols))

        # Electrode particle parameters
        ndD_e[trode]["lambda"] = dD_e[trode]["lambda"]/(k*T_ref)
        ndD_e[trode]["Omga"] = dD_e[trode]["Omga"] / (k*T_ref)
        ndD_e[trode]["Omgb"] = dD_e[trode]["Omgb"] / (k*T_ref)
        ndD_e[trode]["Omgc"] = dD_e[trode]["Omgc"] / (k*T_ref)
        ndD_e[trode]["B"] = dD_e[trode]['B']/(k*T_ref*N_A*dD_e[trode]['cs_ref'])
        ndD_e[trode]["EvdW"] = dD_e[trode]["EvdW"] / (k*T_ref)
        muRfunc = props_am.muRfuncs(ndD_s["T"], ndD_e[trode]).muRfunc
        cs0bar = ndD_s["cs0"][trode]
        cs0 = np.array([cs0bar])
        Type = ndD_e[trode]['type']
        if Type in ndD_s["2varTypes"]:
            cs0 = (cs0, cs0)
            cs0bar = (cs0bar, cs0bar)
        if Type in ndD_s["2varTypes"]:
            ndD_e[trode]["muR_ref"] = -muRfunc(cs0, cs0bar, 0.)[0][0]
        elif Type in ndD_s["1varTypes"]:
            ndD_e[trode]["muR_ref"] = -muRfunc(cs0, cs0bar, 0.)[0]
        ndD_s["phiRef"][trode] = -ndD_e[trode]["muR_ref"][0]
        Nvol, Npart = psd_raw[trode].shape
        # Electrode parameters which depend on the individual
        # particle (e.g. because of non-dimensionalization)
        ndD_e[trode]["indvPart"] = np.empty((Nvol, Npart), dtype=object)
        for i in range(Nvol):
            for j in range(Npart):
                # Bring in a copy of the nondimensional parameters
                # we've calculated so far which are the same for
                # all particles.
                ndD_e[trode]["indvPart"][i,j] = ndD_e[trode].copy()
                # This creates a reference for shorthand.
                ndD_tmp = ndD_e[trode]["indvPart"][i,j]
                # This specific particle dimensions
                ndD_tmp["N"] = ndD_s["psd_num"][trode][i,j]
                plen = dD_s["psd_len"][trode][i,j]
                parea = dD_s["psd_area"][trode][i,j]
                pvol = dD_s["psd_vol"][trode][i,j]
                ndD_tmp["kappa"] = (
                    dD_e[trode]['kappa']
                    / (k*T_ref*N_A*dD_e[trode]['cs_ref']*plen**2))
                ndD_tmp["beta_s"] = (dD_e[trode]['dgammadc']
                                     * plen * N_A*dD_e[trode]['cs_ref']
                                     / dD_e[trode]['kappa'])
                ndD_tmp["D"] = dD_e[trode]['D']*t_ref/plen**2
                ndD_tmp["k0"] = ((parea/pvol)*dD_e[trode]['k0']*t_ref
                                 / (F*dD_e[trode]["cs_ref"]))
                ndD_tmp["Rfilm"] = (
                    dD_e[trode]["Rfilm"] / (
                        t_ref*k*T_ref/(e*F*plen*dD_e[trode]["cs_ref"])))
                ndD_tmp["delta_L"] = pvol/(parea*plen)
                if Type in ["homog_sdn", "homog2_sdn"]:
                    ndD_tmp["Omga"] = size2regsln(plen)

    # Set up macroscopic input information.
    ndDVref = ndD_s["phiRef"]["c"] - ndD_s["phiRef"]["a"]
    # CV setpoint and voltage cutoff values
    ndD_s["Vset"] = -((e/(k*T_ref))*dD_s["Vset"] + ndDVref)
    ndD_s["phimin"] = -((e/(k*T_ref))*dD_s["Vmax"] + ndDVref)
    ndD_s["phimax"] = -((e/(k*T_ref))*dD_s["Vmin"] + ndDVref)

    # Current or voltage segments profiles
    dD_s["segments_tvec"] = np.zeros(2*numsegs + 1)
    dD_s["segments_setvec"] = np.zeros(2*numsegs + 1)
    if ndD_s["profileType"] == "CVsegments":
        dD_s["segments_setvec"][0] = -(k*T_ref/e)*ndDVref
    elif ndD_s["profileType"] == "CCsegments":
        dD_s["segments_setvec"][0] = 0.
    tPrev = 0.
    for segIndx in range(numsegs):
        tNext = tPrev + dD_s["tramp"]
        dD_s["segments_tvec"][2*segIndx+1] = tNext
        tPrev = tNext
        # Factor of 60 here to convert to s
        tNext = tPrev + (segs[segIndx][1] * 60 - dD_s["tramp"])
        dD_s["segments_tvec"][2*segIndx+2] = tNext
        tPrev = tNext
        setNext = segs[segIndx][0]
        dD_s["segments_setvec"][2*segIndx+1] = setNext
        dD_s["segments_setvec"][2*segIndx+2] = setNext
    ndD_s["segments_tvec"] = dD_s["segments_tvec"] / t_ref
    if ndD_s["profileType"] == "CCsegments":
        ndD_s["segments_setvec"] = dD_s["segments_setvec"] / curr_ref
    elif ndD_s["profileType"] == "CVsegments":
        ndD_s["segments_setvec"] = -(
            (e/(k*T_ref))*dD_s["segments_setvec"] + ndDVref)
    if "segments" in ndD_s["profileType"]:
        dD_s["tend"] = dD_s["segments_tvec"][-1]
        # Pad the last segment so no extrapolation occurs
        dD_s["segments_tvec"][-1] = dD_s["tend"]*1.01
    ndD_s["tend"] = dD_s["tend"] / t_ref
    if ndD_s["profileType"] == "CC" and not isClose(ndD_s["currset"], 0.):
        ndD_s["tend"] = np.abs(ndD_s["capFrac"] / ndD_s["currset"])

    return dD_s, ndD_s, dD_e, ndD_e


def distr_part(dD_s, ndD_s, dD_e, ndD_e):
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
            raw = (dD_s["psd_mean"][trode] * np.ones((Nv, Np)))
        else:
            var = stddev**2
            mu = np.log((mean**2)/np.sqrt(var+mean**2))
            sigma = np.sqrt(np.log(var/(mean**2)+1))
            raw = np.random.lognormal(mu, sigma, size=(Nv, Np))
        psd_raw[trode] = raw
        # For particles with internal profiles, convert psd to
        # integers -- number of steps
        solidDisc = dD_e[trode]["disc"]
        if solidType in ["ACR"]:
            psd_num[trode] = (
                np.ceil(psd_raw[trode]/solidDisc).astype(np.integer))
            psd_len[trode] = solidDisc*psd_num[trode]
        elif solidType in ["CHR", "diffn", "CHR2", "diffn2"]:
            psd_num[trode] = (
                np.ceil(psd_raw[trode]/solidDisc).astype(np.integer) + 1)
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
            psd_vol[trode] = (1.2263 * psd_len[trode]**2
                              * dD_e[trode]['thickness'])
        elif solidShape == "cylinder":
            psd_area[trode] = (2 * np.pi * psd_len[trode]
                               * dD_e[trode]['thickness'])
            psd_vol[trode] = (np.pi * psd_len[trode]**2
                              * dD_e[trode]['thickness'])
    return psd_raw, psd_num, psd_len, psd_area, psd_vol


def distr_G(dD, ndD):
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


def size2regsln(size):
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
    if isinstance(param, np.ndarray):
        param[param < 2] = 2.
    else:
        if param < 2:
            param = 2.
    return param


def test_system_input(dD, ndD):
    if not isClose(dD['Tabs'], 298.):
        raise Exception("Temperature dependence not implemented")
    if ndD['Nvol']["c"] < 1:
        raise Exception("Must have at least one porous electrode")
    if ndD["profileType"] not in ["CC", "CV", "CCsegments", "CVsegments"]:
        raise NotImplementedError("profileType {pt} unknown".format(
            pt=ndD["profileType"]))


def test_electrode_input(dD, ndD, dD_s, ndD_s):
    T298 = isClose(dD_s['Tabs'], 298.)
    if not T298:
        raise NotImplementedError("Temperature dependence not supported")
    solidType = ndD['type']
    solidShape = ndD['shape']
    if solidType in ["ACR", "homog_sdn"] and solidShape != "C3":
        raise Exception("ACR and homog_sdn req. C3 shape")
    if (solidType in ["CHR", "diffn"] and solidShape not in
            ["sphere", "cylinder"]):
        raise NotImplementedError("CHR and diffn req. sphere or cylinder")
    if ((solidType not in ndD_s["1varTypes"]) and
            (solidType not in ndD_s["2varTypes"])):
        raise NotImplementedError("Input solidType not defined")
    if solidShape not in ["C3", "sphere", "cylinder"]:
        raise NotImplementedError("Input solidShape not defined")
    if solidType == "homog_sdn" and not T298:
        raise NotImplementedError("homog_snd req. Tabs=298")


def writeConfigFile(P, filename="input_params.cfg"):
    with open(filename, "w") as fo:
        P.write(fo)


def writeDicts(dD, ndD, filenamebase="input_dict"):
    pickle.dump(dD, open(filenamebase + "_dD.p", "wb"))
    pickle.dump(ndD, open(filenamebase + "_ndD.p", "wb"))


def readDicts(filenamebase="input_dict"):
    try:
        dD = pickle.load(open(filenamebase + "_dD.p", "rb"))
        ndD = pickle.load(open(filenamebase + "_ndD.p", "rb"))
    except UnicodeDecodeError:
        dD = pickle.load(open(filenamebase + "_dD.p", "rb"), encoding="latin1")
        ndD = pickle.load(open(filenamebase + "_ndD.p", "rb"), encoding="latin1")
    return dD, ndD


def getConfig(inFile):
    P = configparser.RawConfigParser()
    P.optionxform = str
    P.read(inFile)
    return P


def isClose(a, b):
    if np.abs(a - b) < 1e-12:
        return True
    else:
        return False
