import dash
from dash import html, dcc
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output
import plotly.express as px

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots


from mpet.config import constants

from mpet.plot import plot_data_db

dff, dff_c_sub, dff_cd_sub, dff_bulkp, dff_csld_sub, df_cbar = plot_data_db.main()

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
    if m_select["Config trode type"].iloc[0] in constants.one_var_types:
        trode_type = 1
    elif m_select["Config trode type"].iloc[0] in constants.two_var_types:
        trode_type = 2
    # plots
    subplt_solid_surf_con_c = subplt_solid_surf_con("Cathode", m_select, trode_type)
    subplt_solid_surf_con_a = subplt_solid_surf_con("Anode", m_select, trode_type)
    cbarline_c, cbarline_a = subplots_cbarlinec(m_select, trode_type)
    cbar_c, cbar_c2, cbar_a, cbar_a2 = ani_cbar(select_single_model, trode_type)
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
    range_min = 0.9*min_y if min_y > 0 else 1.1*min_y
    range_max = 1.1*max_y if max_y > 0 else 1.1*max_y
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
                     range=[range_min, range_max])
    fig.update_xaxes(title=u'Battery Position (\u00B5m)')
    fig.update_layout(legend=dict(
        yanchor="bottom",
        y=1.02,
        xanchor="left",
        x=0.01
    ))
    return fig


def subplt_solid_surf_con(trodes, df, trode_type):
    trode = trodes[0].lower()
    try:
        pfx = max(df["pfx"])
        sStr = max(df["sStr"])
        r = int(max(df["Npart"+trode]))
        c = int(max(df["Nvol"+trode]))
        fig = make_subplots(rows=r, cols=c, shared_xaxes=True, shared_yaxes=True,
                            x_title='Time (s)', y_title='Solid surface concentration',
                            row_titles=['Particle ' + str(n) for n in range(1, r+1)],
                            column_titles=['Volume ' + str(n) for n in range(1, c+1)])
        for rr in range(0, r):
            for cc in range(0, c):
                if trode_type == 2:
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
            partStr = "partTrode{trode}vol{vInd}part{pInd}" + max(df["sStr"])
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
                        c1str_base = (max(df["pfx"])
                                      + partStr.format(trode=trode, pInd=rr, vInd=cc) + "c1")
                        c2str_base = (max(df["pfx"])
                                      + partStr.format(trode=trode, pInd=rr, vInd=cc) + "c2")
                        c3str_base = (max(df["pfx"])
                                      + partStr.format(trode=trode, pInd=rr, vInd=cc) + "cav")
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
                        cstr_base = (max(df["pfx"]) + partStr.format(trode=trode, vInd=cc, pInd=rr)
                                     + "c")
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


def subplots_cbarlinec(df, trode_type):
    for trodes in ["Cathode", "Anode"]:
        trode = trodes[0].lower()
        try:
            partStr = ("partTrode{trode}vol{{vInd}}part{{pInd}}".format(trode=trode)
                       + max(df["sStr"]))
            r = int(max(df["Npart"+trode]))
            c = int(max(df["Nvol"+trode]))
            fig = make_subplots(rows=r, cols=c, shared_xaxes=True, shared_yaxes=True,
                                x_title='Time (s)', y_title='Particle Average Filling Fraction',
                                row_titles=['Particle ' + str(n) for n in range(1, r+1)],
                                column_titles=['Volume ' + str(n) for n in range(1, c+1)])
            for rr in range(0, r):
                for cc in range(0, c):
                    datax = df['Time (s)']
                    if trode_type == 2:
                        str1_cbar_base = max(df["pfx"]) + partStr + "c1bar"
                        str2_cbar_base = max(df["pfx"]) + partStr + "c2bar"
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
                        str_cbar_base = max(df["pfx"]) + partStr + "cbar"
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


def ani_cbar(ms_cbar, trode_type):
    def plot(cbar):
        plot = px.scatter(df_cbar_select, x='c', y='r', animation_frame='Time',
                          animation_group='rc',
                          color=cbar,
                          range_color=[0,1],
                          color_continuous_scale='Turbo',
                          size='Relative size',
                          height=None,
                          width=None
                          )
        return plot
    for trodes in ["Cathode", "Anode"]:
        trode = trodes[0].lower()
        df_cbar_select = df_cbar[(df_cbar.Model == ms_cbar) & (df_cbar.Trode == trode)]
        if df_cbar_select.empty:
            cbar = px.line(title='Selected model has no '+trodes.lower())
            cbar2 = px.line()
        else:
            if trode_type == 1:
                cbar = plot('Cbar')
                cbar2 = px.line(title='Selected model has no Cbar2')
            elif trode_type == 2:
                cbar = plot('Cbar1')
                cbar2 = plot('Cbar2')
                cbar2.update_layout(title=trodes, transition_duration=100, height=None, width=None)
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
