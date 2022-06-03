from dash import Dash, dcc, html, Input, Output
import plotly.express as px
import pandas as pd
import os
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

import mpet.geometry as geom
from mpet import mod_cell
from mpet import utils
from mpet.config import Config, constants

#  from plotly.subplots import make_subplots
#  import plotly.graph_objects as go

app = Dash(__name__)

# Define colors to be used
colors = {
    'background': '#FFFFFF',
    'text': '#000000',
    'bg_table': '#B9B6B6',
    'text_table': '#111111',
}

# Import data from all folders in sim_output into one dataframe
# Read in the simulation results and calcuations data
dataDir = "sim_output"
dataFiles = [os.path.join(dataDir, f) for f in os.listdir(dataDir)
             if os.path.isdir(os.path.join(dataDir, f))]
dff = pd.DataFrame()
dff_c_sub = pd.DataFrame()
dff_cd_sub = pd.DataFrame()

for indir in dataFiles:
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
    try:
        ffvec_a = utils.get_dict_key(data, pfx + 'ffrac_a')
        nparta = Npart["a"]
        nvola = Nvol["a"]
    except KeyError:
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
    ylbl_ce = 'Concentration of electrolyte [M]'
    datay_ce = datay_c * c_ref / 1000.
    # elytep
    ylbl_pe = 'Potential of electrolyte [V]'
    datay_pe = datay_p*(k*Tref/e) - Vstd
    cGP_L = utils.get_dict_key(data, "c_lyteGP_L")
    pGP_L = utils.get_dict_key(data, "phi_lyteGP_L")
    cmat = np.hstack((cGP_L.reshape((-1,1)), datay_c, datay_c[:,-1].reshape((-1,1))))
    pmat = np.hstack((pGP_L.reshape((-1,1)), datay_p, datay_p[:,-1].reshape((-1,1))))
    disc = geom.get_elyte_disc(Nvol, config["L"], config["poros"], config["BruggExp"])
    i_edges = np.zeros((numtimes, len(facesvec)))
    for tInd in range(numtimes):
        i_edges[tInd, :] = mod_cell.get_lyte_internal_fluxes(
            cmat[tInd, :], pmat[tInd, :], disc, config)[1]
    # elytei
    ylbl_cd = r'Current density of electrolyte [A/m$^2$]'
    datax_cd = facesvec
    datay_cd = i_edges * (F*constants.c_ref*config["D_ref"]/config["L_ref"])
    # elytedivi
    ylbl_d = r'Divergence of electrolyte current density [A/m$^3$]'
    datax_d = cellsvec
    datay_d = np.diff(i_edges, axis=1) / disc["dxvec"]
    datay_d *= (F*constants.c_ref*config["D_ref"]/config["L_ref"]**2)
    # fraction
    t_current = times
    tfrac = (t_current - tmin)/(tmax - tmin) * 100
    # elytecons
    sep = pfx + 'c_lyte_s'
    anode = pfx + 'c_lyte_a'
    cath = pfx + 'c_lyte_c'
    cvec = utils.get_dict_key(data, cath)
    if Nvol["s"]:
        cvec_s = utils.get_dict_key(data, sep)
        cvec = np.hstack((cvec_s, cvec))
    if "a" in trodes:
        cvec_a = utils.get_dict_key(data, anode)
        cvec = np.hstack((cvec_a, cvec))
    cavg = np.sum(porosvec*dxvec*cvec, axis=1)/np.sum(porosvec*dxvec)

    df = pd.DataFrame({
        "Model": model,
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

    # add data for surface concentration to df
    if "a" in trodes:
        trode = "a"
        str_base = (pfx
                    + "partTrode{trode}vol{{vInd}}part{{pInd}}".format(trode=trode)
                    + sStr + "c")
        for pInd in range(Npart[trode]):
            for vInd in range(Nvol[trode]):
                sol_str = str_base.format(pInd=pInd, vInd=vInd)
                sol_str_data = utils.get_dict_key(data, sol_str, squeeze=False)[:,-1]
                df[sol_str] = sol_str_data
    if "c" in trodes:
        trode = "c"
        str_base = (pfx
                    + "partTrode{trode}vol{{vInd}}part{{pInd}}".format(trode=trode)
                    + sStr + "c")
        for pInd in range(Npart[trode]):
            for vInd in range(Nvol[trode]):
                sol_str = str_base.format(pInd=pInd, vInd=vInd)
                sol_str_data = utils.get_dict_key(data, sol_str, squeeze=False)[:,-1]
                df[sol_str] = sol_str_data

    dff = pd.concat([dff, df], ignore_index=True)

    # build dataframe for plots electrolyte concentration or potential
    dff_c = pd.DataFrame()
    dff_cd = pd.DataFrame()
    for i in range(len(datay_ce)):
        df_c = pd.DataFrame({"Model": model,
                             "fraction": round(tfrac[i]),
                             "fraction orig": tfrac[i],
                             "cellsvec": cellsvec[:],
                             "Concentration electrolyte": datay_ce[i,:],
                             "Potential electrolyte": datay_pe[i,:],
                             "Divergence electrolyte curr dens": datay_d[i,:]
                             })
        dff_c = pd.concat([dff_c, df_c], ignore_index=True)
    for i in range(len(datay_cd)):
        df_cd = pd.DataFrame({"Model": model,
                              "fraction": round(tfrac[i]),
                              "fraction orig": tfrac[i],
                              "facesvec": facesvec[:],
                              "Curreny density electrolyte": datay_cd[i,:]
                              })
        dff_cd = pd.concat([dff_cd, df_cd], ignore_index=True)

    # make subselection dataframe with one fraction per rounded fraction
    for i in np.unique(dff_c["fraction"]):
        df_sub = dff_c[dff_c["fraction"] == i]
        md = 1.0
        mdx = 0.0
        for j in np.unique(df_sub["fraction orig"]):
            dx = abs(i-j)
            if dx < md:
                md = dx
                mdx = j
        select = dff_c[dff_c["fraction orig"] == mdx]
        dff_c_sub = pd.concat([dff_c_sub, select], ignore_index=True)
############################
# Define plots
############################
# fig = px.line(dff, x="Time (s)", y="Cathode Filling Fraction", color="Model")

# fig2 = px.line(dff, x="Time (s)", y="Anode Filling Fraction", color="Model")

# Define Markdown text
markdown_text = '''
# MPET visualisation

This dashboard shows visualisations of all the MPET simulation output saved in the folder
"sim_output."
'''

defaultmodel = (dff['Model'].unique()[0])
# Define components of app
app.layout = html.Div(style={'backgroundColor': colors['background']},
                      children=[
                      html.Div(style={'backgroundColor': 'DodgerBlue'},
                               children=[
                               dcc.Markdown(
                                   children=markdown_text,
                                   style={'textAlign': 'center',
                                          'color': colors['text'],
                                          'font-family':'Sans-serif'}),
                                   html.H4(children='Select models to display in all plots:',
                                           style={'font-family':'Sans-serif'}),
                                   dcc.Checklist(options=dff['Model'].unique(),
                                                 value=dff['Model'].unique(),
                                                 id='model-selection',
                                                 labelStyle={'display': 'block'},
                                                 style={'font-family':'Sans-serif',
                                                        'margin-bottom': '20px',
                                                        'background': 'DodgerBlue'}),
                                   html.Hr(style={"color": 'black', 'borderWidth': '10'})]),

                      html.H3(children='Voltage',
                              style={'textAlign': 'center', 'font-family':'Sans-serif'}),
                      dcc.Dropdown(['Time (s)', 'Cathode Filling Fraction'], 'Time (s)',
                                   id='xaxis-column',
                                   style={'width':'50%', 'font-family':'Sans-serif'}),
                      dcc.Graph(id='Voltage-graph-double'),

                      html.Hr(style={"color": 'black'}),
                      html.Hr(style={"color": 'black'}),

                      html.H3(children='Solid surface concentration anode',
                              style={'textAlign': 'center', 'font-family':'Sans-serif'}),
                      dcc.Dropdown(options=dff['Model'].unique(), value=defaultmodel,
                                   id='select_surf_a_model',
                                   style={'width':'50%', 'font-family':'Sans-serif'}),
                      dcc.Graph(id='Surface-concentration-anode'),

                      html.H3(children='Solid surface concentration cathode',
                              style={'textAlign': 'center', 'font-family':'Sans-serif'}),
                      dcc.Dropdown(options=dff['Model'].unique(), value=defaultmodel,
                                   id='select_surf_c_model',
                                   style={'width':'50%', 'font-family':'Sans-serif'}),
                      dcc.Graph(id='Surface-concentration-cathode'),

                      html.Hr(style={"color": 'black'}),
                      html.Hr(style={"color": 'black'}),

                      html.H3(children='Overall utilization / state of charge of electrode',
                              style={'textAlign': 'center', 'font-family':'Sans-serif'}),
                      dcc.Graph(id='Cathode-filling-fraction'),
                      dcc.Graph(id='Anode-filling-fraction'),

                      html.Hr(style={"color": 'black'}),
                      html.Hr(style={"color": 'black'}),

                      html.H3(children='Average concentration of electrolyte',
                              style={'textAlign': 'center', 'font-family':'Sans-serif'}),
                      dcc.Graph(id='elytecons'),

                      html.Hr(style={"color": 'black'}),
                      html.Hr(style={"color": 'black'}),

                      html.H3(children='Electrolyte concentration',
                              style={'textAlign': 'center', 'font-family':'Sans-serif'}),

                      html.Hr(style={"color": 'black'}),
                      html.Hr(style={"color": 'black'}),

                      html.H3(children='Current profile',
                              style={'textAlign': 'center', 'font-family':'Sans-serif'}),
                      dcc.Graph(id="current"),

                      html.Hr(style={"color": 'black'}),
                      html.Hr(style={"color": 'black'}),

                      html.H3(children='Power',
                              style={'textAlign': 'center', 'font-family':'Sans-serif'}),
                      dcc.Graph(id="power"),

                      html.Hr(style={"color": 'black'}),
                      html.Hr(style={"color": 'black'}),

                      html.H3(children='Electrolyte concentration and potential',
                              style={'textAlign': 'center', 'font-family':'Sans-serif'}),
                      dcc.Graph(id='electrolyte-concentration-ani'),

                      html.Hr(style={"color": 'black'}),
                      html.Hr(style={"color": 'black'}),

                      html.H3(children='Solid particle-average concentrations',
                              style={'textAlign': 'center', 'font-family':'Sans-serif'}),

                      html.Hr(style={"color": 'black'}),
                      html.Hr(style={"color": 'black'}),

                      html.H3(children='All solid concentrations or potentials',
                              style={'textAlign': 'center', 'font-family':'Sans-serif'}),

                      html.Hr(style={"color": 'black'}),
                      html.Hr(style={"color": 'black'}),

                      html.H3(children='Average solid concentrations',
                              style={'textAlign': 'center', 'font-family':'Sans-serif'}),

                      html.Hr(style={"color": 'black'}),
                      html.Hr(style={"color": 'black'}),

                      html.H3(children='Cathode potential',
                              style={'textAlign': 'center', 'font-family':'Sans-serif'}),

                      html.Hr(style={"color": 'black'}),
                      html.Hr(style={"color": 'black'})
                      ])


# Do alllll the callbacks
@app.callback(
    Output('Voltage-graph-double', 'figure'),
    Input('xaxis-column', 'value'),
    Input('model-selection', 'value')
    )
def update_graph(xaxis_column_name, model_selection
                 ):

    fig = px.line(x=dff[xaxis_column_name][np.in1d(dff['Model'], model_selection)],
                  y=dff['Voltage (V)'][np.in1d(dff['Model'], model_selection)],
                  color=dff["Model"][np.in1d(dff['Model'], model_selection)])
    fig.update_layout(margin={'l': 40, 'b': 80, 't': 10, 'r': 0}, hovermode='closest')
    fig.update_yaxes(title='Voltage (V)',
                     type='linear')
    fig.update_xaxes(title=xaxis_column_name,
                     type='linear'
                     )
    return fig


@app.callback(
    Output('Surface-concentration-anode', 'figure'),
    Input('select_surf_a_model', 'value')
    )
def update_graph_surfa(select_model):
    str_base = (pfx
                + "partTrode{trode}vol{{vInd}}part{{pInd}}".format(trode='a')
                + sStr + "c")
    r = int(max(dff["Nparta"][dff['Model'] == select_model]))
    c = int(max(dff["Nvola"][dff['Model'] == select_model]))
    fig = make_subplots(rows=r, cols=c, shared_xaxes=True, shared_yaxes=True,
                        x_title='Time [s]', y_title='Solid surface concentration',
                        row_titles=['Particle ' + str(n) for n in range(1, r+1)],
                        column_titles=['Volume ' + str(n) for n in range(1, c+1)])
    for rr in range(0, r):
        for cc in range(0, c):
            sol_str = str_base.format(pInd=rr, vInd=cc)
            datax = dff['Time (s)'][dff['Model'] == select_model]
            datay = dff[sol_str][dff['Model'] == select_model]
            fig.add_trace(
                trace=go.Scatter(x=datax, y=datay, line_color='darkslategray'),
                row=rr+1, col=cc+1)
    fig.update_yaxes(range=[0,1.01])
    fig.update_layout(showlegend=False)
    return fig


@app.callback(
    Output('Surface-concentration-cathode', 'figure'),
    Input('select_surf_c_model', 'value')
    )
def update_graph_surfc(select_model):
    str_base = (pfx
                + "partTrode{trode}vol{{vInd}}part{{pInd}}".format(trode='c')
                + sStr + "c")
    r = int(max(dff["Npartc"][dff['Model'] == select_model]))
    c = int(max(dff["Nvolc"][dff['Model'] == select_model]))
    fig = make_subplots(rows=r, cols=c, shared_xaxes=True, shared_yaxes=True,
                        x_title='Time (s)', y_title='Solid surface concentration',
                        row_titles=['Particle ' + str(n) for n in range(1, r+1)],
                        column_titles=['Volume ' + str(n) for n in range(1, c+1)])
    for rr in range(0, r):
        for cc in range(0, c):
            sol_str = str_base.format(pInd=rr, vInd=cc)
            datax = dff['Time (s)'][dff['Model'] == select_model]
            datay = dff[sol_str][dff['Model'] == select_model]
            fig.add_trace(
                trace=go.Scatter(x=datax, y=datay, line_color='darkslategray'),
                row=rr+1, col=cc+1)
    fig.update_yaxes(range=[0,1.01])
    fig.update_layout(showlegend=False)
    return fig


@app.callback(
    Output('Cathode-filling-fraction', 'figure'),
    Input('model-selection', 'value')
    )
def update_graph2(model_selection):
    fig = px.line(x=dff['Time (s)'][np.in1d(dff['Model'], model_selection)],
                  y=dff['Cathode Filling Fraction'][np.in1d(dff['Model'], model_selection)],
                  color=dff["Model"][np.in1d(dff['Model'], model_selection)])
    fig.update_yaxes(range=[0,1], title='Cathode Filling Fraction')
    fig.update_xaxes(title='Time [s]')

    fig.update_layout(margin={'l': 40, 'b': 80, 't': 10, 'r': 0}, hovermode='closest')

    return fig


@app.callback(
    Output('Anode-filling-fraction', 'figure'),
    Input('model-selection', 'value')
    )
def update_graph3(model_selection):

    fig = px.line(x=dff['Time (s)'][np.in1d(dff['Model'], model_selection)],
                  y=dff['Anode Filling Fraction'][np.in1d(dff['Model'], model_selection)],
                  color=dff["Model"][np.in1d(dff['Model'], model_selection)])
    fig.update_yaxes(range=[0,1], title='Anode Filling Fraction')
    fig.update_xaxes(title='Time [s]')
    fig.update_layout(margin={'l': 40, 'b': 80, 't': 10, 'r': 0}, hovermode='closest')
    return fig


@app.callback(
    Output('elytecons', 'figure'),
    Input('model-selection', 'value')
    )
def update_graph32(model_selection):
    fig = px.line(x=dff['Time (s)'][np.in1d(dff['Model'], model_selection)],
                  y=dff['cavg'][np.in1d(dff['Model'], model_selection)],
                  color=dff["Model"][np.in1d(dff['Model'], model_selection)])
    fig.update_yaxes(title='Avg. Concentration of electrolyte [nondim]')
    fig.update_xaxes(title='Time [s]')
    fig.update_layout(margin={'l': 40, 'b': 80, 't': 10, 'r': 0}, hovermode='closest')
    return fig


@app.callback(
    Output('current', 'figure'),
    Input('model-selection', 'value')
    )
def update_graph4(model_selection):

    fig = px.line(x=dff['Time (s)'][np.in1d(dff['Model'], model_selection)],
                  y=dff['Current'][np.in1d(dff['Model'], model_selection)],
                  color=dff["Model"][np.in1d(dff['Model'], model_selection)])
    fig.update_yaxes(title='Current [C-rate]')
    fig.update_xaxes(title='Time [s]')
    fig.update_layout(margin={'l': 40, 'b': 80, 't': 10, 'r': 0}, hovermode='closest')
    return fig


@app.callback(
    Output('power', 'figure'),
    Input('model-selection', 'value')
    )
def update_graph5(model_selection):

    fig = px.line(x=dff['Time (s)'][np.in1d(dff['Model'], model_selection)],
                  y=dff['Power'][np.in1d(dff['Model'], model_selection)],
                  color=dff["Model"][np.in1d(dff['Model'], model_selection)])
    fig.update_yaxes(title=u'Power [W/m\u00b2]')
    fig.update_xaxes(title='Time [s]')
    fig.update_layout(margin={'l': 40, 'b': 80, 't': 10, 'r': 0}, hovermode='closest')
    return fig


@app.callback(
    Output('electrolyte-concentration-ani', 'figure'),
    Input('model-selection', 'value'))
def display_animated_graph(model_selection):
    #  print('model_selection',model_selection)
    #  print(dff_c[round(dff_c['fraction'], 1) == selected_perc])
    #  filtered_df = dff_c[np.in1d(dff['Model'], model_selection)]
    #  filtered_df = dff_c[round(dff_c['fraction'], 1) == selected_perc]
    m_select = dff_c_sub[np.in1d(dff_c_sub['Model'], model_selection)]
    max_y = max(m_select["Concentration electrolyte"])
    min_y = min(m_select["Concentration electrolyte"])
    fig = px.line(dff_c_sub[np.in1d(dff_c_sub['Model'], model_selection)],
                  x="cellsvec",
                  y="Concentration electrolyte",
                  color="Model",
                  animation_frame="fraction")
    fig.update_yaxes(title='Concentration of electrolyte [M]',
                     range=[0.9*min_y, 1.1*max_y])
    fig.update_xaxes(title=' Battery Position [microm] ')
    fig.update_layout(margin={'l': 40, 'b': 80, 't': 10, 'r': 0}, hovermode='closest')
    fig.update_layout(transition_duration=500)

    return fig


if __name__ == '__main__':
    app.run_server(debug=True)
