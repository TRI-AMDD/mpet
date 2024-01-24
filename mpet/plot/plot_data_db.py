import pandas as pd
import numpy as np
import os
import argparse
from argparse import RawTextHelpFormatter

import mpet.geometry as geom
from mpet import mod_cell
from mpet import utils
from mpet.config import Config, constants
from mpet.exceptions import UnknownParameterError


###################################################################
# Part 1
# Import data from all subfolders in provided dataDir
# Prepare all data (get everything in correct dataframe) for plotting
# This is based on plot_data.py
###################################################################
def main():
    desc = """ Dashboard that shows all plots and compares the resutls of different models."""
    parser = argparse.ArgumentParser(description=desc, formatter_class=RawTextHelpFormatter)
    parser.add_argument('-d', '--dataDir',
                        help='Directory that contains subfolders with simulation output')
    args = parser.parse_args()

    dataFiles = [os.path.join(args.dataDir, f) for f in os.listdir(args.dataDir)
                 if os.path.isdir(os.path.join(args.dataDir, f))]
    dff = pd.DataFrame()
    dff_c_sub = pd.DataFrame()
    dff_cd_sub = pd.DataFrame()
    dff_csld_sub = pd.DataFrame()
    dff_bulkp = pd.DataFrame()
    df_cbar = pd.DataFrame()

    for indir in dataFiles:
        print('Loading data from:', indir)
        pfx = 'mpet.'
        sStr = "_"
        # Read in the simulation results and calcuations data
        model = os.path.basename(indir)
        dataFileName = "output_data"
        dataFile = os.path.join(indir, dataFileName)
        data = utils.open_data_file(dataFile)
        try:
            utils.get_dict_key(data, pfx + 'current')
        except KeyError:
            pfx = ''
        try:
            utils.get_dict_key(data, pfx + "partTrodecvol0part0" + sStr + "cbar")
        except KeyError:
            sStr = "."
        # Read in the parameters used to define the simulation
        config = Config.from_dicts(indir)
        # simulated (porous) electrodes
        trodes = config["trodes"]
        # Pick out some useful calculated values
        limtrode = config["limtrode"]
        k = constants.k                      # Boltzmann constant, J/(K Li)
        Tref = constants.T_ref               # Temp, K
        e = constants.e                      # Charge of proton, C
        F = constants.F                      # C/mol
        c_ref = constants.c_ref
        td = config["t_ref"]
        Etheta = {"a": 0.}
        cap = config[limtrode, "cap"]
        for trode in trodes:
            Etheta[trode] = -(k*Tref/e) * config[trode, "phiRef"]
        Vstd = Etheta["c"] - Etheta["a"]
        Nvol = config["Nvol"]
        Npart = config["Npart"]
        psd_len = config["psd_len"]
        # Discretization (and associated porosity)
        Lfac = 1e6
        dxc = config["L"]["c"]/Nvol["c"]
        dxvec = np.array(Nvol["c"] * [dxc])
        porosvec = np.array(Nvol["c"] * [config["poros"]["c"]])
        cellsvec = dxc*np.arange(Nvol["c"]) + dxc/2.
        if config["have_separator"]:
            dxs = config["L"]["s"]/Nvol["s"]
            dxvec_s = np.array(Nvol["s"] * [dxs])
            dxvec = np.hstack((dxvec_s, dxvec))
            poros_s = np.array(Nvol["s"] * [config["poros"]["s"]])
            porosvec = np.hstack((poros_s, porosvec))
            cellsvec += config["L"]["s"] / config["L"]["c"]
            cellsvec_s = dxs*np.arange(Nvol["s"]) + dxs/2.
            cellsvec = np.hstack((cellsvec_s, cellsvec))
        if "a" in trodes:
            dxa = config["L"]["a"]/Nvol["a"]
            dxvec_a = np.array(Nvol["a"] * [dxa])
            dxvec = np.hstack((dxvec_a, dxvec))
            poros_a = np.array(Nvol["a"] * [config["poros"]["a"]])
            porosvec = np.hstack((poros_a, porosvec))
            cellsvec += config["L"]["a"] / config["L"]["c"]
            cellsvec_a = dxa*np.arange(Nvol["a"]) + dxa/2.
            cellsvec = np.hstack((cellsvec_a, cellsvec))
        cellsvec *= config["L_ref"] * Lfac
        facesvec = np.insert(np.cumsum(dxvec), 0, 0.) * config["L_ref"] * Lfac
        # Extract the reported simulation times
        times = utils.get_dict_key(data, pfx + 'phi_applied_times')
        numtimes = len(times)
        tmin = np.min(times)
        tmax = np.max(times)
        # Voltage profile
        timestd = times*td
        voltage = (Vstd - (k*Tref/e)*utils.get_dict_key(data, pfx + 'phi_applied'))
        # surface concentration
        # soc profile
        ffvec_c = utils.get_dict_key(data, pfx + 'ffrac_c')
        if "a" in trodes:
            ffvec_a = utils.get_dict_key(data, pfx + 'ffrac_a')
            nparta = Npart["a"]
            nvola = Nvol["a"]
        else:
            ffvec_a = 0
            nparta = 0
            nvola = 0
        # Elytecons
        # current
        theoretical_1C_current = config[config['limtrode'], "cap"] / 3600.  # A/m^2
        current = (utils.get_dict_key(data, pfx + 'current')
                   * theoretical_1C_current / config['1C_current_density'] * config['curr_ref'])
        # Power
        current_p = utils.get_dict_key(data, pfx + 'current') * (3600/td) * (cap/3600)  # in A/m^2
        voltage_p = (Vstd - (k*Tref/e)*utils.get_dict_key(data, pfx + 'phi_applied'))  # in V
        power = np.multiply(current_p, voltage_p)

        # Electrolyte concetration / potential
        datax = cellsvec
        c_sep, p_sep = pfx + 'c_lyte_s', pfx + 'phi_lyte_s'
        c_anode, p_anode = pfx + 'c_lyte_a', pfx + 'phi_lyte_a'
        c_cath, p_cath = pfx + 'c_lyte_c', pfx + 'phi_lyte_c'
        datay_c = utils.get_dict_key(data, c_cath, squeeze=False)
        datay_p = utils.get_dict_key(data, p_cath, squeeze=False)
        L_c = config['L']["c"] * config['L_ref'] * Lfac
        Ltot = L_c
        if config["have_separator"]:
            datay_s_c = utils.get_dict_key(data, c_sep, squeeze=False)
            datay_s_p = utils.get_dict_key(data, p_sep, squeeze=False)
            datay_c = np.hstack((datay_s_c, datay_c))
            datay_p = np.hstack((datay_s_p, datay_p))
            L_s = config['L']["s"] * config['L_ref'] * Lfac
            Ltot += L_s
        else:
            L_s = 0
        if "a" in trodes:
            datay_a_c = utils.get_dict_key(data, c_anode, squeeze=False)
            datay_a_p = utils.get_dict_key(data, p_anode, squeeze=False)
            datay_c = np.hstack((datay_a_c, datay_c))
            datay_p = np.hstack((datay_a_p, datay_p))
            L_a = config['L']["a"] * config['L_ref'] * Lfac
            Ltot += L_a
        else:
            L_a = 0
        # elytec
        datay_ce = datay_c * c_ref / 1000.
        # elytep
        datay_pe = datay_p*(k*Tref/e) - Vstd
        i_edges = np.zeros((numtimes, len(facesvec)))
        # elytei & elytedivi
        try:
            cGP_L = utils.get_dict_key(data, "c_lyteGP_L")
            pGP_L = utils.get_dict_key(data, "phi_lyteGP_L")
            cmat = np.hstack((cGP_L.reshape((-1,1)), datay_c, datay_c[:,-1].reshape((-1,1))))
            pmat = np.hstack((pGP_L.reshape((-1,1)), datay_p, datay_p[:,-1].reshape((-1,1))))
            disc = geom.get_elyte_disc(Nvol, config["L"], config["poros"], config["BruggExp"])
            for tInd in range(numtimes):
                i_edges[tInd, :] = mod_cell.get_lyte_internal_fluxes(
                    cmat[tInd, :], pmat[tInd, :], disc, config)[1]
            datay_cd = i_edges * (F*constants.c_ref*config["D_ref"]/config["L_ref"])
            datay_d = np.diff(i_edges, axis=1) / disc["dxvec"]
            datay_d *= (F*constants.c_ref*config["D_ref"]/config["L_ref"]**2)
        except (UnknownParameterError, KeyError):
            datay_cd = i_edges
            datay_d = np.zeros((numtimes, len(cellsvec)))
            datay_d *= (F*constants.c_ref*config["D_ref"]/config["L_ref"]**2)
        # fraction
        t_current = times
        tfrac = (t_current - tmin)/(tmax - tmin) * 100
        # elytecons
        sep = pfx + 'c_lyte_s'
        anode = pfx + 'c_lyte_a'
        cath = pfx + 'c_lyte_c'
        cvec = utils.get_dict_key(data, cath)
        if config["have_separator"]:
            cvec_s = utils.get_dict_key(data, sep)
            cvec = np.hstack((cvec_s, cvec))
        if "a" in trodes:
            cvec_a = utils.get_dict_key(data, anode)
            cvec = np.hstack((cvec_a, cvec))
        try:
            cavg = np.sum(porosvec*dxvec*cvec, axis=1)/np.sum(porosvec*dxvec)
        except np.AxisError:
            cavg = np.sum(porosvec*dxvec*cvec, axis=0)/np.sum(porosvec*dxvec)
        # Get all data in dataframes
        df = pd.DataFrame({
            "Model": model,
            "sStr": sStr,
            "pfx": pfx,
            "Config trode type": config[trode, "type"],
            "Voltage (V)": voltage,
            "Cathode Filling Fraction": ffvec_c,
            "Anode Filling Fraction": ffvec_a,
            "Time (s)": timestd,
            "Npartc": Npart["c"],
            "Nvolc": Nvol["c"],
            "Nparta": nparta,
            "Nvola": nvola,
            "Current": current,
            "Power": power,
            "cavg": cavg
        })
        for trode in trodes:
            partStr = "partTrode{trode}vol{{vInd}}part{{pInd}}".format(trode=trode) + sStr
            lens = psd_len[trode]
            size_fracs = 0.4*np.ones((Nvol[trode], Npart[trode]))
            if np.max(lens) != np.min(lens):
                size_fracs = (lens - np.min(lens))/(np.max(lens) - np.min(lens))
            size_frac_min = 0.2
            sizes = (size_fracs*(1-size_frac_min) + size_frac_min) / Nvol[trode]
            for pInd in range(Npart[trode]):
                for vInd in range(Nvol[trode]):
                    # for c data and cbar data
                    if config[trode, "type"] in constants.one_var_types:
                        str_base = (pfx + partStr + "c")
                        sol_str = str_base.format(pInd=pInd, vInd=vInd)
                        sol_str_data = utils.get_dict_key(data, sol_str, squeeze=False)[:,-1]
                        str_cbar_base = pfx + partStr + "cbar"
                        sol_cbar_str = str_cbar_base.format(pInd=pInd, vInd=vInd)
                        sol_cbar_str_data = utils.get_dict_key(data, sol_cbar_str)
                        df = pd.concat((df, pd.DataFrame({sol_str: sol_str_data,
                                                          sol_cbar_str: sol_cbar_str_data})),
                                       axis=1)
                        # for cbar movie
                        df_c = pd.DataFrame({
                            'Model': model,
                            "Config trode type": 1,
                            'Time': np.round(timestd),
                            'Cbar': sol_cbar_str_data,
                            'r': pInd,
                            'c': vInd,
                            'rc': str(pInd)+str(vInd),
                            'Relative size': sizes[vInd,pInd],
                            'Trode': trode
                            })
                    elif config[trode, "type"] in constants.two_var_types:
                        str1_base = (pfx + partStr + "c1")
                        str2_base = (pfx + partStr + "c2")
                        sol1_str = str1_base.format(pInd=pInd, vInd=vInd)
                        sol2_str = str2_base.format(pInd=pInd, vInd=vInd)
                        sol1_str_data = utils.get_dict_key(data, sol1_str, squeeze=False)[:,-1]
                        sol2_str_data = utils.get_dict_key(data, sol2_str, squeeze=False)[:,-1]
                        str1_cbar_base = pfx + partStr + "c1bar"
                        str2_cbar_base = pfx + partStr + "c2bar"
                        sol1_cbar_str = str1_cbar_base.format(pInd=pInd, vInd=vInd)
                        sol2_cbar_str = str2_cbar_base.format(pInd=pInd, vInd=vInd)
                        sol1_cbar_str_data = utils.get_dict_key(data, sol1_cbar_str)
                        sol2_cbar_str_data = utils.get_dict_key(data, sol2_cbar_str)
                        df = pd.concat((df, pd.DataFrame({sol1_str: sol1_str_data,
                                                          sol2_str: sol2_str_data,
                                                          sol1_cbar_str: sol1_cbar_str_data,
                                                          sol2_cbar_str: sol2_cbar_str_data})),
                                       axis=1)
                        df_c = pd.DataFrame({
                            'Model': model,
                            "Config trode type": 2,
                            'Time': np.round(timestd),
                            'Cbar1': sol1_cbar_str_data,
                            'Cbar2': sol2_cbar_str_data,
                            'r': pInd,
                            'c': vInd,
                            'rc': str(pInd)+str(vInd),
                            'Relative size': sizes[vInd,pInd],
                            'Trode': trode
                            })
                    df_cbar = pd.concat([df_cbar, df_c])
        dff = pd.concat([dff, df], ignore_index=True)
        # build dataframe for plots electrolyte concentration or potential
        # and for csld subplot animation (time, pind, vind, y)
        dff_c = pd.DataFrame({"Model": model,
                              "Time fraction": np.round(np.repeat(tfrac, np.shape(datay_ce)[1])),
                              "fraction orig": np.repeat(tfrac, np.shape(datay_ce)[1]),
                              "cellsvec": np.tile(cellsvec, np.shape(datay_ce)[0]),
                              "Concentration electrolyte": datay_ce.flatten(),
                              "Potential electrolyte": datay_pe.flatten(),
                              "Divergence electrolyte curr dens": datay_d.flatten()
                              })
        dff_cd = pd.DataFrame({"Model": model,
                               "Time fraction": np.round(np.repeat(tfrac, np.shape(datay_cd)[1])),
                               "fraction orig": np.repeat(tfrac, np.shape(datay_cd)[1]),
                               "facesvec": np.tile(facesvec, np.shape(datay_cd)[0]),
                               "Curreny density electrolyte": datay_cd.flatten()
                               })
        # Build dataframes for bulkp and csld
        dff_bulkp_c = pd.DataFrame()
        dff_csld = pd.DataFrame()
        # cstr can have varying length, determine maximum length
        if config[trode, "type"] in constants.one_var_types:
            partStr = "partTrode{trode}vol{vInd}part{pInd}" + sStr
            cstr_base = pfx + partStr + "c"
            maxlength = max([np.shape(utils.get_dict_key(data, cstr_base.format(
                            trode=t, pInd=p, vInd=v)))[1]
                for t in trodes for p in range(Npart[t]) for v in range(Nvol[t])])
        else:
            partStr = "partTrode{trode}vol{vInd}part{pInd}" + sStr
            maxlength = max([np.shape(utils.get_dict_key(data, (pfx + partStr + "c1").format(
                            trode=t, pInd=p, vInd=v)))[1]
                for t in trodes for p in range(Npart[t]) for v in range(Nvol[t])])
        for trode in trodes:
            bulkp = pfx + 'phi_bulk_{trode}'.format(trode=trode)
            dataybulkp = utils.get_dict_key(data, bulkp).flatten()
            if trode == "a":
                dataxbulkp = cellsvec[:Nvol["a"]]
            elif trode == "c":
                dataxbulkp = cellsvec[-Nvol["c"]:]
            datat = np.repeat(timestd, len(dataxbulkp))
            datatfrac = np.repeat(tfrac, len(dataxbulkp))
            dataxbulkp = np.repeat([dataxbulkp], len(timestd), axis=0).flatten()
            df_b = pd.DataFrame({
                "Model": model,
                "Time (s)": datat,
                "Time fraction (%)": np.round(datatfrac),
                "fraction orig": datatfrac,
                "Trode": trode,
                "Potential (nondim)": dataybulkp,
                "Position in electrode": dataxbulkp})
            dff_bulkp_c = pd.concat([dff_bulkp_c, df_b])
            partStr = "partTrode{trode}vol{{vInd}}part{{pInd}}".format(trode=trode) + sStr
            for pInd in range(Npart[trode]):
                for vInd in range(Nvol[trode]):
                    lens_str = "lens_{vInd}_{pInd}".format(vInd=vInd, pInd=pInd)
                    if config[trode, "type"] in constants.one_var_types:
                        cstr_base = pfx + partStr + "c"
                        cstr = cstr_base.format(trode=trode, pInd=pInd, vInd=vInd)
                        datay = np.empty([len(timestd), maxlength]) + np.nan
                        yy = utils.get_dict_key(data, cstr)
                        datay[0:len(timestd), 0:np.shape(yy)[1]] = yy
                        datax = np.empty(maxlength) + np.nan
                        datax[0:np.shape(yy)[1]] = np.linspace(0, psd_len[trode][vInd,pInd] * Lfac,
                                                               np.shape(yy)[1])
                        if trode == trodes[0] and pInd == 0 and vInd == 0:
                            df_csld = pd.DataFrame({"Model": model,
                                                    "sStr": sStr,
                                                    "pfx": pfx,
                                                    "Config trode type": config[trode, "type"],
                                                    "Npartc": Npart["c"],
                                                    "Nvolc": Nvol["c"],
                                                    "Nparta": nparta,
                                                    "Nvola": nvola,
                                                    "time (s)": np.repeat(timestd,
                                                                          np.shape(datay)[1]),
                                                    "Time fraction": np.repeat(np.round(tfrac),
                                                                               np.shape(datay)[1]),
                                                    "fraction orig": np.repeat(tfrac,
                                                                               np.shape(datay)[1]),
                                                    lens_str: np.repeat([datax], len(datay),
                                                                        axis=0).flatten(),
                                                    cstr: datay.flatten()
                                                    })
                        else:
                            if lens_str not in df_csld.columns.to_numpy():
                                lens_str_data = np.repeat([datax], len(datay), axis=0).flatten()
                                df_csld = pd.concat((df_csld, pd.DataFrame({lens_str:
                                                                           lens_str_data})),
                                                    axis=1)
                            df_csld = pd.concat((df_csld, pd.DataFrame({cstr: datay.flatten()})),
                                                axis=1)
                    elif config[trode, "type"] in constants.two_var_types:
                        c1str_base = pfx + partStr + "c1"
                        c2str_base = pfx + partStr + "c2"
                        c3str_base = pfx + partStr + "cav"
                        c1str = c1str_base.format(trode=trode, pInd=pInd, vInd=vInd)
                        c2str = c2str_base.format(trode=trode, pInd=pInd, vInd=vInd)
                        c3str = c3str_base.format(trode=trode, pInd=pInd, vInd=vInd)
                        datay1 = datay2 = datay3 = np.empty([len(timestd), maxlength]) + np.nan
                        yy1 = utils.get_dict_key(data, c1str)
                        datay1[0:len(timestd), 0:np.shape(yy1)[1]] = yy1
                        yy2 = utils.get_dict_key(data, c2str)
                        datay2[0:len(timestd), 0:np.shape(yy2)[1]] = yy2
                        datay3 = 0.5*(datay1 + datay2)
                        datax = np.empty(maxlength) + np.nan
                        numy = np.shape(yy1)[1] if isinstance(yy1, np.ndarray) else 1
                        datax[0:np.shape(yy1)[1]] = np.linspace(0,
                                                                psd_len[trode][vInd,pInd] * Lfac,
                                                                numy)
                        if trode == trodes[0] and pInd == 0 and vInd == 0:
                            df_csld = pd.DataFrame({"Model": model,
                                                    "sStr": sStr,
                                                    "pfx": pfx,
                                                    "Config trode type": config[trode, "type"],
                                                    "Npartc": Npart["c"],
                                                    "Nvolc": Nvol["c"],
                                                    "Nparta": nparta,
                                                    "Nvola": nvola,
                                                    "time (s)": np.repeat(timestd,
                                                                          np.shape(datay1)[1]),
                                                    "Time fraction":
                                                        np.repeat(np.round(tfrac),
                                                                  np.shape(datay1)[1]),
                                                    "fraction orig":
                                                        np.repeat(tfrac,
                                                                  np.shape(datay1)[1]),
                                                    lens_str: np.repeat([datax], len(datay1),
                                                                        axis=0).flatten(),
                                                    c1str: datay1.flatten(),
                                                    c2str: datay2.flatten(),
                                                    c3str: datay3.flatten()
                                                    })
                        else:
                            if lens_str not in df_csld.columns.to_numpy():
                                df_csld[lens_str] = np.repeat([datax], len(datay1),
                                                              axis=0).flatten(),
                            df_csld[c1str] = datay1.flatten()
                            df_csld[c2str] = datay2.flatten()
                            df_csld[c3str] = datay3.flatten()
            dff_csld = pd.concat([dff_csld, df_csld], ignore_index=True)
        # make subselection dataframe with one fraction per rounded fraction
        for i in np.unique(dff_c["Time fraction"]):
            df_sub = dff_c[dff_c["Time fraction"] == i]
            md = 1.0
            mdx = 0.0
            for j in np.unique(df_sub["fraction orig"]):
                dx = abs(i-j)
                if dx < md:
                    md = dx
                    mdx = j
            select = dff_c[dff_c["fraction orig"] == mdx]
            dff_c_sub = pd.concat([dff_c_sub, select], ignore_index=True)
            select = dff_cd[dff_cd["fraction orig"] == mdx]
            dff_cd_sub = pd.concat([dff_cd_sub, select], ignore_index=True)
            select = dff_csld[dff_csld["fraction orig"] == mdx]
            dff_csld_sub = pd.concat([dff_csld_sub, select], ignore_index=True)
            select = dff_bulkp_c[dff_bulkp_c["fraction orig"] == mdx]
            dff_bulkp = pd.concat([dff_bulkp, select], ignore_index=True)
    return dff, dff_c_sub, dff_cd_sub, dff_bulkp, dff_csld_sub, df_cbar
