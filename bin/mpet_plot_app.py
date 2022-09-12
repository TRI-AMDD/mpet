import dash
from dash import html, dcc
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output
import plotly.express as px
import pandas as pd
import os
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
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
    ttl_fmt = "% = {perc:2.1f}"
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
    dataReporter = config["dataReporter"]
    Nvol = config["Nvol"]
    Npart = config["Npart"]
    psd_len = config["psd_len"]
    # Discretization (and associated porosity)
    Lfac = 1e6
    Lunit = r"$\mu$m"
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
    xmin = 0
    xmax = Ltot
    # elytec
    datay_ce = datay_c * c_ref / 1000.
    # elytep
    datay_pe = datay_p*(k*Tref/e) - Vstd
    i_edges = np.zeros((numtimes, len(facesvec)))
    datax_cd = facesvec
    datax_d = cellsvec
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
        spacing = 1.0 / Nvol[trode]
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
            for t in trodes for p in range(Npart[trode]) for v in range(Nvol[trode])])
    else:
        partStr = "partTrode{trode}vol{vInd}part{pInd}" + sStr
        maxlength = max([np.shape(utils.get_dict_key(data, (pfx + partStr + "c1").format(
                        trode=t, pInd=p, vInd=v)))[1]
            for t in trodes for p in range(Npart[trode]) for v in range(Nvol[trode])])
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
                                                "Config trode type": config[trode, "type"],
                                                "Npartc": Npart["c"],
                                                "Nvolc": Nvol["c"],
                                                "Nparta": nparta,
                                                "Nvola": nvola,
                                                "time (s)": np.repeat(timestd, np.shape(datay)[1]),
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
                            df_csld = pd.concat((df_csld, pd.DataFrame({lens_str: lens_str_data})),
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
                    datax[0:np.shape(yy1)[1]] = np.linspace(0, psd_len[trode][vInd,pInd] * Lfac,
                                                            numy)
                    if trode == trodes[0] and pInd == 0 and vInd == 0:
                        df_csld = pd.DataFrame({"Model": model,
                                                "Config trode type": config[trode, "type"],
                                                "Npartc": Npart["c"],
                                                "Nvolc": Nvol["c"],
                                                "Nparta": nparta,
                                                "Nvola": nvola,
                                                "time (s)": np.repeat(timestd,
                                                                      np.shape(datay1)[1]),
                                                "Time fraction": np.repeat(np.round(tfrac),
                                                                           np.shape(datay1)[1]),
                                                "fraction orig": np.repeat(tfrac,
                                                                           np.shape(datay1)[1]),
                                                lens_str: np.repeat([datax], len(datay1),
                                                                    axis=0).flatten(),
                                                c1str: datay1.flatten(),
                                                c2str: datay2.flatten(),
                                                c3str: datay3.flatten()
                                                })
                    else:
                        if lens_str not in df_csld.columns.to_numpy():
                            df_csld[lens_str] = np.repeat([datax], len(datay1), axis=0).flatten(),
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

###################################################################
# Part 2
# Define plots (define that there is a grpah and perform callback (see part 3))
# style elements
###################################################################
app = dash.Dash(__name__, external_stylesheets=[dbc.themes.CERULEAN])

# Define colors to be used
colors = {
    'background': '#FFFFFF',
    'lithium': '#2fa4e7',
    'text': '#000000',
    'bg_table': '#B9B6B6',
    'text_table': '#111111',
}
# Define Markdown text
markdown_text = '''
# MPET visualisation

This dashboard shows visualisations of all the MPET simulation output saved in the folder
''' + args.dataDir

defaultmodel = (dff['Model'].unique()[0])
# Define components of app
app.layout = html.Div([
    dbc.Row(dbc.Col(
        html.H1("MPET visualisation",
            style={'textAlign': 'center', 'background': colors['lithium'], 'color': 'white',
                   'margin-top': '10px', 'margin-bottom': '10px'}),
        ),),
    dbc.Row(dbc.Col([
        html.Div("Select models to display in all plots",
            style={'color': colors['lithium'], 'margin-left': '20px'}),
        dcc.Checklist(
            options=dff['Model'].unique(),
            value=dff['Model'].unique(),
            id='model-selection',
            labelStyle={'display': 'block'},
            style={'margin-bottom': '20px', 'margin-left': '20px',},
            )]
        ),),
    html.Hr(style={"color": colors['lithium'], "height":'5px'}),
    html.H2("General properties", style={'margin-left': '20px'}),
    html.Hr(style={"color": colors['lithium'], "height":'5px'}),
    dbc.Row([
        dbc.Col([
            # v/vt
            dbc.Card([
                html.H4(
                    children='Voltage',
                    style={'textAlign': 'center', 'font-family':'Sans-serif', 'margin-top': '10px'}
                    ),
                dcc.Dropdown(
                    ['Time (s)', 'Cathode Filling Fraction'],
                    'Time (s)',
                    id='xaxis-column',
                    style={'font-family':'Sans-serif', 'margin-left': '10px', 'width': '80%'}),
                dcc.Loading(children=[dcc.Graph(id='Voltage-graph-double')],
                            color=colors['lithium'], type="dot", fullscreen=False)],
                style={'margin-bottom': '20px', 'margin-left': '10px'}
                ),
            # power
            dbc.Card([
                html.H4(children='Power', style={
                    'textAlign': 'center', 'font-family':'Sans-serif', 'margin-top': '10px'}),
                dcc.Loading(children=[dcc.Graph(id="power")],
                            color=colors['lithium'], type="dot", fullscreen=False,)
                ], style={'margin-bottom': '20px', 'margin-left': '10px'})
            ],
            width=4,),
        dbc.Col([
            # Curr
            dbc.Card([
                html.H4(
                    children='Current profile',
                    style={'textAlign': 'center', 'font-family':'Sans-serif',
                           'margin-bottom': '45px', 'margin-top': '10px'}),
                dcc.Loading(children=[dcc.Graph(id="current")],
                            color=colors['lithium'], type="dot", fullscreen=False),
            ], style={'margin-bottom': '20px'}),
            # elytecons
            dbc.Card([
                html.H4(children='Average concentration of electrolyte',
                        style={'textAlign': 'center', 'font-family':'Sans-serif',
                               'margin-top': '10px'}),
                dcc.Loading(children=[dcc.Graph(id='elytecons'),],
                            color=colors['lithium'], type="dot", fullscreen=False,),
            ], style={'margin-bottom': '20px'})],
            width=4,),
        dbc.Col([
            # soc
            dbc.Card([
                html.H4(children='Overall utilization',
                        style={'textAlign': 'center', 'font-family':'Sans-serif',
                               'margin-top': '10px'}),
                html.H4(children='(state of charge of electrode)',
                        style={'textAlign': 'center', 'font-family':'Sans-serif'}),
                dcc.Loading(children=[
                    html.H5(
                        children='Cathode',
                        style={'textAlign': 'center', 'font-family':'Sans-serif',
                            'color': '#7f7f7f', 'margin-top': '10px'}),
                    dcc.Graph(id='Cathode-filling-fraction'),
                    html.H5(
                        children='Anode',
                        style={'textAlign': 'center', 'font-family':'Sans-serif',
                               'color': '#7f7f7f', 'margin-top': '5px', 'margin-bottom': '0px'}),
                    dcc.Graph(id='Anode-filling-fraction'),
                    ],
                    color=colors['lithium'], type="dot", fullscreen=False,),
                ], style={'margin-bottom': '20px', 'margin-right': '10px'})],
                width=4,),
        ]
    ),
    html.Hr(style={"color": colors['lithium'], "height":'5px'}),
    html.H2("Details of electrolyte", style={'margin-left': '20px'}),
    html.Hr(style={"color": colors['lithium'], "height":'5px'}),
    dbc.Row([
        dbc.Col([
            # elytec
            dbc.Card([
                html.H4(children='Electrolyte concentration',
                        style={'textAlign': 'center', 'margin-top': '10px'}),
                dcc.Loading(children=[dcc.Graph(id='electrolyte-concentration-ani')],
                            color=colors['lithium'], type="dot", fullscreen=False,),
            ], style={'margin-bottom': '20px', 'margin-left': '10px'}),
            # elytep
            dbc.Card([
                html.H4(children='Electrolyte potential',
                        style={'textAlign': 'center', 'margin-top': '10px'}),
                dcc.Loading(children=[dcc.Graph(id='electrolyte-potential-ani'),],
                            color=colors['lithium'], type="dot", fullscreen=False,),
            ], style={'margin-bottom': '20px', 'margin-left': '10px'}),
        ], width=6),
        dbc.Col([
            # elytei
            dbc.Card([
                html.H4(children='Electrolyte current density',
                        style={'textAlign': 'center', 'margin-top': '10px'}),
                dcc.Loading(children=[dcc.Graph(id='electrolyte-cd-ani'),],
                            color=colors['lithium'], type="dot", fullscreen=False,),
            ], style={'margin-bottom': '20px', 'margin-right': '10px'}),
            # elytedivi
            dbc.Card([
                html.H4(children='Divergence of electrolyte current density',
                        style={'textAlign': 'center', 'margin-top': '10px'}),
                dcc.Loading(children=[dcc.Graph(id='electrolyte-decd-ani'),],
                            color=colors['lithium'], type="dot", fullscreen=False,),
            ], style={'margin-bottom': '20px', 'margin-right': '10px'}),
        ], width=6),
    ]),
    # bulkp
    dbc.Row(dbc.Col(dbc.Card([
        html.H4(children='Macroscopic cathode/anode solid phase potential',
                style={'textAlign': 'center', 'margin-top': '10px'}),
        dcc.Loading(children=[dcc.Graph(id='bulkp'),],
                    color=colors['lithium'], type="dot", fullscreen=False,),
        ], style={'margin-left': '10px', 'margin-right': '10px'}), width=12)),
    #
    html.Hr(style={"color": colors['lithium'], "height":'5px'}),
    html.H2("Details of concentrations for individual models", style={'margin-left': '20px'}),
    html.Hr(style={"color": colors['lithium'], "height":'5px'}),
    dbc.Row(dbc.Col([
        html.Div("Select model to display",
            style={'color': colors['lithium'], 'margin-left': '20px'}),
        dcc.Dropdown(options=dff['Model'].unique(), value=defaultmodel,
                     id='select_single_model',
                     style={'width':'50%', 'margin-left': '10px', 'margin-bottom': '20px'})]
        ),),
    # cbar
    dbc.Row(dbc.Col(dbc.Card([
        html.H4(children='Average solid concentrations',
                style={'textAlign': 'center', 'font-family':'Sans-serif', 'margin-top': '10px'}),
        dcc.Loading(children=[dcc.Graph(id='cbar_c'),
                              dcc.Graph(id='cbar_c2'),
                              dcc.Graph(id='cbar_a'),
                              dcc.Graph(id='cbar_a2'),],
                    color=colors['lithium'], type="dot", fullscreen=False,),
    ], style={'margin-left': '10px', 'margin-right': '10px', 'margin-bottom': '20px'}))),
    # surf
    dbc.Row(dbc.Col(dbc.Card([
        html.H4(children='Solid surface concentration',
                style={'textAlign': 'center', 'font-family':'Sans-serif', 'margin-top': '10px'}),
        dcc.Loading(children=[dcc.Graph(id='Surface-concentration-cathode'),
                              dcc.Graph(id='Surface-concentration-anode'),],
                    color=colors['lithium'], type="dot", fullscreen=False,),
    ], style={'margin-left': '10px', 'margin-right': '10px', 'margin-bottom': '20px',
              "overflow": "scroll"}))),
    # csld
    dbc.Row(dbc.Col(dbc.Card([
        html.H4(children='All solid concentrations',
                style={'textAlign': 'center', 'font-family':'Sans-serif', 'margin-top': '10px'}),
        html.H5(children='Time percentage',
                style={'textAlign': 'left', 'font-family':'Sans-serif', 'margin-left': '10px'}),
        dcc.Slider(
            0,100,step=1,value=0,
            id='timefraction_slider',
            marks={
                0: '0', 5: '5', 10: '10', 15: '15', 20: '20',
                25: '25', 30: '30', 35: '35', 40: '40', 45: '45',
                50: '50', 55: '55', 60: '60', 65: '65', 70: '70',
                75: '75', 80: '80', 85: '85', 90: '90', 95: '95',
                100: '100'
            }, tooltip={"placement": "top", "always_visible": True}),
        dcc.Graph(id='csld_c'),
        dcc.Graph(id='csld_a'),
    ], style={'margin-left': '10px', 'margin-right': '10px', 'margin-bottom': '20px',
              "overflow": "scroll"}))),
    # cbarline
    dbc.Row(dbc.Col(dbc.Card([
        html.H4(children='Average concentration in each particle of electrode',
                style={'textAlign': 'center', 'font-family':'Sans-serif', 'margin-top': '10px'}),
        dcc.Loading(children=[dcc.Graph(id='cbarline_c'),
                              dcc.Graph(id='cbarline_a'),],
                    color=colors['lithium'], type="dot", fullscreen=False,),
    ], style={'margin-left': '10px', 'margin-right': '10px', 'margin-bottom': '20px',
              "overflow": "scroll"}))),
])


###################################################################
# Part 3
# Define callback functions to get interactive plots
###################################################################
@app.callback(
    Output('current', 'figure'),
    Output('power', 'figure'),
    Output('Cathode-filling-fraction', 'figure'),
    Output('Anode-filling-fraction', 'figure'),
    Output('elytecons', 'figure'),
    Input('model-selection', 'value')
    )
def update_graphs_multimodels(model_selection):
    m_select = dff[np.in1d(dff['Model'], model_selection)]
    # plots
    current = plot_multimodel(m_select, 'Current', ytitle='Current (C-rate)')
    power = plot_multimodel(m_select, 'Power', ytitle=u'Power (W/m\u00b2)')
    fillfrac_c = plot_multimodel(m_select, 'Cathode Filling Fraction',
                                 ytitle='Cathode Filling Fraction')
    fillfrac_a = plot_multimodel(m_select, 'Anode Filling Fraction',
                                 ytitle='Anode Filling Fraction')
    elytecons = plot_multimodel(m_select, 'cavg',
                                ytitle='Avg. Concentration of electrolyte (nondim)')
    return current, power, fillfrac_c, fillfrac_a, elytecons


@app.callback(
    Output('Voltage-graph-double', 'figure'),
    Input('xaxis-column', 'value'),
    Input('model-selection', 'value')
    )
def update_graphs_multimodels_voltage(xaxis_column_name, model_selection):
    m_select = dff[np.in1d(dff['Model'], model_selection)]
    voltage = plot_multimodel(m_select, 'Voltage (V)', xaxis=xaxis_column_name)
    return voltage


@app.callback(
    Output('electrolyte-concentration-ani', 'figure'),
    Output('electrolyte-potential-ani', 'figure'),
    Output('electrolyte-cd-ani', 'figure'),
    Output('electrolyte-decd-ani', 'figure'),
    Output('bulkp', 'figure'),
    Input('model-selection', 'value')
    )
def callback_multimodel_movies(model_selection):
    m_select_dff_c = dff_c_sub[np.in1d(dff_c_sub['Model'], model_selection)]
    m_select_dff_cd = dff_cd_sub[np.in1d(dff_cd_sub['Model'], model_selection)]
    m_select_dff_bulkp = dff_bulkp[np.in1d(dff_bulkp['Model'], model_selection)]
    electrolyte_concentration = ani_elytrolyte(m_select_dff_c, "cellsvec",
                                               "Concentration electrolyte",
                                               "Time fraction", 'Concentration of electrolyte (M)')
    eletrolyte_potential = ani_elytrolyte(m_select_dff_c, "cellsvec", "Potential electrolyte",
                                          "Time fraction", 'Potential of electrolyte (V)')
    curr_dens = ani_elytrolyte(m_select_dff_cd, "facesvec",
                               "Curreny density electrolyte",
                               "Time fraction", u"Current density of electrolyte (A/m\u00b2)")
    div_curr_dens = ani_elytrolyte(m_select_dff_c, "cellsvec",
                                   "Divergence electrolyte curr dens",
                                   "Time fraction",
                                   u"Divergence electrolyte current density (A/m\u00b3)")
    bulkp = ani_bulkp(m_select_dff_bulkp)
    return electrolyte_concentration, eletrolyte_potential, curr_dens, div_curr_dens, bulkp


@app.callback(
    Output('Surface-concentration-cathode', 'figure'),
    Output('Surface-concentration-anode', 'figure'),
    Output('cbarline_c', 'figure'),
    Output('cbarline_a', 'figure'),
    Output('cbar_c', 'figure'),
    Output('cbar_c2', 'figure'),
    Output('cbar_a', 'figure'),
    Output('cbar_a2', 'figure'),
    Input('select_single_model', 'value')
    )
def update_graphs_single_models(select_single_model):
    m_select = dff[dff['Model'] == select_single_model]
    # plots
    subplt_solid_surf_con_c = subplt_solid_surf_con("Cathode", m_select)
    subplt_solid_surf_con_a = subplt_solid_surf_con("Anode", m_select)
    cbarline_c, cbarline_a = subplots_cbarlinec(m_select)
    cbar_c, cbar_c2, cbar_a, cbar_a2 = ani_cbar(select_single_model)
    return (subplt_solid_surf_con_c, subplt_solid_surf_con_a,
            cbarline_c, cbarline_a, cbar_c, cbar_c2, cbar_a, cbar_a2
            )


@app.callback(
    Output('csld_c', 'figure'),
    Output('csld_a', 'figure'),
    Input('select_single_model', 'value'),
    Input('timefraction_slider', 'value')
    )
def update_csld(select_single_model, timefraction_slider):
    m_select_csld = dff_csld_sub[dff_csld_sub['Model'] == select_single_model]
    csld_c, csld_a = subplots_t_csld(m_select_csld, timefraction_slider)
    return (csld_c, csld_a)


# Plot functions
def plot_multimodel(df, yaxis, ytitle=None, xaxis='Time (s)'):
    plot = px.line(df, x=xaxis, y=yaxis, color="Model")
    if xaxis != 'Time (s)':
        plot.update_xaxes(title=xaxis)
    if ytitle:
        plot.update_yaxes(title=ytitle)
    plot.update_layout(legend=dict(
        yanchor="bottom",
        y=1.02,
        xanchor="left",
        x=0.01
    ))
    return plot


def ani_elytrolyte(df, xname, yname, ani, ytitle):
    max_y = max(df[yname])
    min_y = min(df[yname])
    if not (max_y == 0 and min_y == 0):
        fig = px.line(df,
                      x=xname,
                      y=yname,
                      color="Model",
                      animation_frame=ani,
                      animation_group=xname)
    else:
        fig = px.line(title='Data not availale for selected model(s)')
    fig.update_yaxes(title=ytitle,
                     range=[0.9*min_y, 1.1*max_y])
    fig.update_xaxes(title=u'Battery Position (\u00B5m)')
    fig.update_layout(legend=dict(
        yanchor="bottom",
        y=1.02,
        xanchor="left",
        x=0.01
    ))
    return fig


def subplt_solid_surf_con(trodes, df):
    trode = trodes[0].lower()
    try:
        if df["Config trode type"].iloc[0] in constants.one_var_types:
            type2c = False
        elif df["Config trode type"].iloc[0] in constants.two_var_types:
            type2c = True
        r = int(max(df["Npart"+trode]))
        c = int(max(df["Nvol"+trode]))
        fig = make_subplots(rows=r, cols=c, shared_xaxes=True, shared_yaxes=True,
                            x_title='Time (s)', y_title='Solid surface concentration',
                            row_titles=['Particle ' + str(n) for n in range(1, r+1)],
                            column_titles=['Volume ' + str(n) for n in range(1, c+1)])
        for rr in range(0, r):
            for cc in range(0, c):
                if type2c:
                    str1_base = (pfx
                                 + "partTrode{trode}vol{{vInd}}part{{pInd}}".format(trode=trode)
                                 + sStr + "c1")
                    str2_base = (pfx
                                 + "partTrode{trode}vol{{vInd}}part{{pInd}}".format(trode=trode)
                                 + sStr + "c2")
                    sol1_str = str1_base.format(pInd=rr, vInd=cc)
                    sol2_str = str2_base.format(pInd=rr, vInd=cc)
                    datax = df['Time (s)']
                    datay1 = df[sol1_str]
                    datay2 = df[sol2_str]
                    fig.add_trace(
                        trace=go.Scatter(x=datax, y=datay1, line_color='red', name='c1'),
                        row=rr+1, col=cc+1)
                    fig.add_trace(
                        trace=go.Scatter(x=datax, y=datay2, line_color='blue', name='c2'),
                        row=rr+1, col=cc+1)
                else:
                    str_base = (pfx
                                + "partTrode{trode}vol{{vInd}}part{{pInd}}".format(trode=trode)
                                + sStr + "c")
                    sol_str = str_base.format(pInd=rr, vInd=cc)
                    datax = df['Time (s)']
                    datay = df[sol_str]
                    fig.add_trace(
                        trace=go.Scatter(x=datax, y=datay, line_color='darkslategray'),
                        row=rr+1, col=cc+1)
        fig.update_yaxes(range=[0,1.01])
        fig.update_layout(height=((r+1)*150), width=((c+1)*150), showlegend=False, title=trodes)
    except ValueError:
        fig = px.line(title='Selected model has no '+trodes.lower())
        fig.update_xaxes(visible=False, showgrid=False)
        fig.update_yaxes(visible=False, showgrid=False)
    return fig


def subplots_t_csld(df, tf):
    for trodes in ["Cathode", "Anode"]:
        try:
            trode = trodes[0].lower()
            partStr = "partTrode{trode}vol{vInd}part{pInd}" + sStr
            r = int(max(df["Npart"+trode]))
            c = int(max(df["Nvol"+trode]))
            fig = make_subplots(rows=r, cols=c, shared_xaxes=True, shared_yaxes=True,
                                x_title=u'Position (\u00B5m)',
                                y_title='Solid Concentrations of Particles in Electrode',
                                row_titles=['Particle ' + str(n) for n in range(1, r+1)],
                                column_titles=['Volume ' + str(n) for n in range(1, c+1)])
            type2c = False
            df_select = df[df['Time fraction'] == int(tf)]
            if df_select["Config trode type"].iloc[0] in constants.one_var_types:
                type2c = False
            elif df_select["Config trode type"].iloc[0] in constants.two_var_types:
                type2c = True
            for rr in range(0, r):
                for cc in range(0, c):
                    lens_str = "lens_{vInd}_{pInd}".format(vInd=cc, pInd=rr)
                    if type2c is True:
                        c1str_base = pfx + partStr.format(trode=trode, pInd=rr, vInd=cc) + "c1"
                        c2str_base = pfx + partStr.format(trode=trode, pInd=rr, vInd=cc) + "c2"
                        c3str_base = pfx + partStr.format(trode=trode, pInd=rr, vInd=cc) + "cav"
                        datax = df_select[lens_str]
                        datay1 = df_select[c1str_base]
                        datay2 = df_select[c2str_base]
                        datay3 = df_select[c3str_base]
                        fig.add_trace(
                            trace=go.Scatter(x=datax, y=datay1, line_color='red', name='c1'),
                            row=rr+1, col=cc+1)
                        fig.add_trace(
                            trace=go.Scatter(x=datax, y=datay2, line_color='blue', name='c2'),
                            row=rr+1, col=cc+1)
                        fig.add_trace(
                            trace=go.Scatter(x=datax, y=datay3, line_color='grey',
                                             name='avg. c1 & c2'),
                            row=rr+1, col=cc+1)
                    else:
                        cstr_base = pfx + partStr.format(trode=trode, vInd=cc, pInd=rr) + "c"
                        cstr = cstr_base
                        datax = df_select[lens_str]
                        datay = df_select[cstr]
                        fig.add_trace(
                            trace=go.Scatter(x=datax, y=datay, line_color='darkslategray'),
                            row=rr+1, col=cc+1)
            fig.update_yaxes(range=[0,1.01])
            fig.update_layout(height=((r+1)*150), width=((c+1)*150), showlegend=False,
                              title=trodes+', time = {time} s'.format(
                              time=round(df_select["time (s)"].iloc[0])))
        except ValueError:
            fig = px.line(title='Selected model has no '+trodes.lower())
            fig.update_xaxes(visible=False, showgrid=False)
            fig.update_yaxes(visible=False, showgrid=False)
        if trode == "c":
            fig1 = fig
    return fig1, fig


def subplots_cbarlinec(df):
    for trodes in ["Cathode", "Anode"]:
        trode = trodes[0].lower()
        try:
            partStr = "partTrode{trode}vol{{vInd}}part{{pInd}}".format(trode=trode) + sStr
            r = int(max(df["Npart"+trode]))
            c = int(max(df["Nvol"+trode]))
            fig = make_subplots(rows=r, cols=c, shared_xaxes=True, shared_yaxes=True,
                                x_title='Time (s)', y_title='Particle Average Filling Fraction',
                                row_titles=['Particle ' + str(n) for n in range(1, r+1)],
                                column_titles=['Volume ' + str(n) for n in range(1, c+1)])
            # this does not work if models with multiple plot types
            if (df["Config trode type"].iloc[0] in constants.one_var_types):
                type2c = False
                str_cbar_base = pfx + partStr + "cbar"
            elif (df["Config trode type"].iloc[0] in constants.two_var_types):
                type2c = True
                str1_cbar_base = pfx + partStr + "c1bar"
                str2_cbar_base = pfx + partStr + "c2bar"
            for rr in range(0, r):
                for cc in range(0, c):
                    datax = df['Time (s)']
                    if type2c is True:
                        sol1_str = str1_cbar_base.format(pInd=rr, vInd=cc)
                        sol2_str = str2_cbar_base.format(pInd=rr, vInd=cc)
                        datay = df[sol1_str]
                        datay2 = df[sol2_str]
                        fig.add_trace(
                            trace=go.Scatter(x=datax, y=datay, line_color='red', name='c1bar'),
                            row=rr+1, col=cc+1)
                        fig.add_trace(
                            trace=go.Scatter(x=datax, y=datay2, line_color='blue', name='c2bar'),
                            row=rr+1, col=cc+1)
                    else:
                        sol_str = str_cbar_base.format(pInd=rr, vInd=cc)
                        datay = df[sol_str]
                        fig.add_trace(
                            trace=go.Scatter(x=datax, y=datay, line_color='darkslategray'),
                            row=rr+1, col=cc+1)
            fig.update_yaxes(range=[0,1.01])
            fig.update_layout(height=((r+1)*150), width=((c+1)*150),
                              showlegend=False, title=trodes)
        except ValueError:
            fig = px.line(title='Selected model has no '+trodes.lower())
            fig.update_xaxes(visible=False, showgrid=False)
            fig.update_yaxes(visible=False, showgrid=False)
        if trode == "c":
            fig1 = fig
    return fig1, fig


def ani_cbar(ms_cbar):
    for trodes in ["Cathode", "Anode"]:
        trode = trodes[0].lower()
        df_cbar_select = df_cbar[(df_cbar.Model == ms_cbar) & (df_cbar.Trode == trode)]
        if df_cbar_select.empty:
            cbar = cbar2 = px.line(title='Selected model has no '+trodes.lower())
            cbar.update_xaxes(visible=False, showgrid=False)
            cbar.update_yaxes(visible=False, showgrid=False)
        else:
            if df_cbar_select["Config trode type"].iloc[0] == 1:
                cbar = px.scatter(df_cbar_select, x='c', y='r', animation_frame='Time',
                                  animation_group='rc',
                                  color='Cbar',
                                  range_color=[0,1],
                                  color_continuous_scale='Turbo',
                                  size='Relative size'
                                  )
                cbar2 = px.scatter()
            elif df_cbar_select["Config trode type"].iloc[0] == 2:
                cbar = px.scatter(df_cbar_select, x='c', y='r', animation_frame='Time',
                                  animation_group='rc',
                                  color='Cbar1',
                                  range_color=[0,1],
                                  color_continuous_scale='Turbo',
                                  size='Relative size'
                                  )
                cbar2 = px.scatter(df_cbar_select, x='c', y='r', animation_frame='Time',
                                   animation_group='rc',
                                   color='Cbar2',
                                   range_color=[0,1],
                                   color_continuous_scale='Turbo',
                                   size='Relative size'
                                   )
                cbar2.update_layout(title=trodes, transition_duration=100)
            cbar.update_layout(title=trodes, transition_duration=100)
        cbar.update_xaxes(visible=False, showgrid=False)
        cbar.update_yaxes(visible=False, showgrid=False)
        cbar2.update_xaxes(visible=False, showgrid=False)
        cbar2.update_yaxes(visible=False, showgrid=False)
        if trode == "c":
            cbar_c = cbar
            cbar_c2 = cbar2
    return cbar_c, cbar_c2, cbar, cbar2


def ani_bulkp(df):
    bulkp = px.line(df, x='Position in electrode', y='Potential (nondim)',
                    animation_frame='Time fraction (%)', color='Model',
                    line_dash='Trode', log_y=True)
    bulkp.update_yaxes(range=[-4, 1.1*np.log10(max(df['Potential (nondim)']))])
    bulkp.update_xaxes(title='Position in electrode (\u00B5m)')
    return bulkp


if __name__ == '__main__':
    app.run_server(debug=True)
