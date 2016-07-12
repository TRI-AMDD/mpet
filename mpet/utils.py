import numpy as np

import daetools.pyDAE as dae


def set_plot_defaults(mpl):
    axtickfsize = 18
    labelfsize = 20
    legfsize = labelfsize - 2
    txtfsize = labelfsize - 2
    lwidth = 3
    markersize = 10
    markeredgewidth = 0.1
    mpl.rcParams['xtick.labelsize'] = axtickfsize
    mpl.rcParams['ytick.labelsize'] = axtickfsize
    mpl.rcParams['axes.labelsize'] = labelfsize
    mpl.rcParams['font.size'] = txtfsize
    mpl.rcParams['legend.fontsize'] = legfsize
    mpl.rcParams['lines.linewidth'] = lwidth
    mpl.rcParams['lines.markersize'] = markersize
    mpl.rcParams['lines.markeredgewidth'] = markeredgewidth
#    mpl.rcParams['text.usetex'] = True


def mean_linear(a, b=None):
    if isinstance(a, np.ndarray):
        return 0.5*(a[1:] + a[:-1])
    else:
        return 0.5*(a + b)


def mean_harmonic(a, b=None):
    if isinstance(a, np.ndarray):
        return (2 * a[1:] * a[:-1]) / (a[1:] + a[:-1] + 1e-20)
    else:
        return (2 * a * b) / (a + b + 1e-20)


def get_cell_Ntot(Nvol):
    """Nvol is a dictionary containing the number of volumes in each simulated battery section."""
    return np.sum(list(Nvol.values()))


def add_gp_to_vec(vec):
    """Add ghost points to the beginning and end of a vector for applying boundary conditions."""
    out = np.empty(len(vec) + 2, dtype=object)
    out[1:-1] = vec
    return out


def pad_vec(vec):
    """Repeat a vector's first and last values, extending its length by two."""
    out = add_gp_to_vec(vec)
    out[0] = out[1]
    out[-1] = out[-2]
    return out


def get_const_vec(val, N):
    out = np.array([val for indx in range(N)], dtype=object)
    return out


def get_var_vec(var, N, dt=False):
    if dt is True:
#        out = np.array([var.dt(indx) for indx in range(var.NumberOfPoints)])
        out = np.array([var.dt(indx) for indx in range(N)])
    else:
#        out = np.array([var(indx) for indx in range(var.NumberOfPoints)])
        out = np.array([var(indx) for indx in range(N)])
    return out


#def get_asc_const_vec(var, Nvol):
#    """Get a numpy array for a parameter spanning the anode, separator, and cathode."""
##    out = np.zeros(get_cell_Ntot(Nvol))
##    out[0:Nvol["a"]] = [var["a"] for vInd in range(Nvol["a"])]  # anode
##    out[Nvol["a"]:Nvol["a"]+Nvol["s"]] = [var["s"] for vInd in range(Nvol["s"])]  # separator
##    out[Nvol["a"]+Nvol["s"]:] = [var["c"] for vInd in range(Nvol["c"])]  # cathode
#    varout = {}
#    for sectn in ["a", "s", "c"]:
#        # If we have information within this battery section
#        if sectn in var.keys():
#            varout[sectn] = get_const_vec(var[sectn], Nvol[sectn])
#        # Otherwise, fill with zeros
#        else:
#            varout[sectn] = np.zeros(Nvol[sectn])
#    out = np.hstack((varout["a"], varout["s"], varout["c"]))
##    out = np.hstack((get_const_vec(var["a"], Nvol["a"]),
##                     get_const_vec(var["s"], Nvol["s"]),
##                     get_const_vec(var["c"], Nvol["c"])))
#    return out


#def get_asc_var_vec(var, Nvol, dt=False):
def get_asc_vec(var, Nvol, dt=False):
    """Get a numpy array for a variable spanning the anode, separator, and cathode."""
    varout = {}
    for sectn in ["a", "s", "c"]:
        # If we have information within this battery section
        if sectn in var.keys():
            # If it's an array of dae variable objects
            if isinstance(var[sectn], dae.pyCore.daeVariable):
                varout[sectn] = get_var_vec(var[sectn], Nvol[sectn], dt)
            # Otherwise, it's a parameter that varies with electrode section
            else:
                varout[sectn] = get_const_vec(var[sectn], Nvol[sectn])
        # Otherwise, fill with zeros
        else:
            varout[sectn] = np.zeros(Nvol[sectn])
    out = np.hstack((varout["a"], varout["s"], varout["c"]))
#    out = np.empty(get_cell_Ntot(Nvol), dtype=object)
#    if Nvol["a"] > 0 and "a" in var.keys():
#        var_a = get_var_vec(var["a"], dt)
#    else:
#        var_a = np.array([])
#    if Nvol["s"] > 0 and "s" in var.keys():
#        var_s = get_var_vec(var["s"], dt)
#    else:
#        var_s = np.array([])
#    var_c = get_var_vec(var["s"], dt)
#    out = np.hstack((var_a, var_s, var_c))
#    out = np.hstack((get_var_vec(var["a"], dt),
#                     get_var_vec(var["s"], dt),
#                     get_var_vec(var["c"], dt)))
    return out


#def get_elyte_varvec(var, Nvol, dt=False):
#    Nlyte = np.sum(list(Nvol.values()))
#    out = np.empty(Nlyte, dtype=object)
#    # Anode
#    if dt is False:
#        out[0:Nvol["a"]] = [var["a"](vInd) for vInd in range(Nvol["a"])]
#    else:
#        out[0:Nvol["a"]] = [var["a"].dt(vInd) for vInd in range(Nvol["a"])]
#    # Separator: If not present, fill with zeros
#    if Nvol["s"] and "s" in var.keys():
#        if dt is False:
#            out[Nvol["a"]:Nvol["a"]+Nvol["s"]] = [var["s"](vInd) for vInd in range(Nvol["s"])]
#        else:
#            out[Nvol["a"]:Nvol["a"]+Nvol["s"]] = [var["s"].dt(vInd) for vInd in range(Nvol["s"])]
#    else:
#        out[Nvol["a"]:Nvol["a"] + Nvol["s"]] = [0. for vInd in range(Nvol["s"])]
#    # Cathode
#    if dt is False:
#        out[Nvol["a"] + Nvol["s"]:Nlyte] = [var["c"](vInd) for vInd in range(Nvol["c"])]
#    else:
#        out[Nvol["a"] + Nvol["s"]:Nlyte] = [var["c"].dt(vInd) for vInd in range(Nvol["c"])]
#    return out


#def get_padded_asc_vec(var, Nvol):
#    """Get a numpy array for a variable whose value only depends on region (anode, separator, and
#    cathode) with "padding" on the ends (repeating the first and final values).
#    """
#    out = get_asc_vec_with_gp(var, Nvol)
#    out[0] = out[1]
#    out[-1] = out[-2]
##    out = np.empty(np.sum(list(Nvol.values())) + 2, dtype=object)
##    out[0:Nvol["a"]+1] = [var["a"] for vInd in range(Nvol["a"]+1)]
##    out[Nvol["a"]+1:Nvol["a"]+1 + Nvol["s"]] = [var["s"] for vInd in range(Nvol["s"])]
##    out[Nvol["a"]+1 + Nvol["s"]:] = [var["c"] for vInd in range(Nvol["c"]+1)]
#    return out


#def get_asc_vec_with_gp(var, Nvol):
#    """Get a numpy array with ghost points for a variable spanning the anode, separator, and
#    cathode.
#    """
#    tmp = get_asc_vec(var, Nvol)
#    out = np.empty(len(tmp) + 2, dtype=object)
#    out[1:-1] = tmp
#    return out


#def get_padded_vec(var, N=None):
#    """Get a numpy vector of a variable or float padded with its end member values in the ghost
#    points. If var is a float instead of an array, N specifies the output length.
#    """
#    out = get_vec_with_gp(var, N)
#    out[0] = out[1]
#    out[-1] = out[-2]
#    return out


#def get_vec_with_gp(var, N=None):
#    """Get a numpy array with empty ghost points at the edges from a float, a dae variable array,
#    or another numpy array. N should be passed only if var is a float.
#    """
#    if N is None:
#        N = len(var)
#    else:
#        out = np.empty(N + 2)
#        out[1:-1] = [var for indx in range(N)]
#    if isinstance(var, dae.pyCore.daeVariable):
#        out = np.empty(N + 2, dtype=object)
#        out[1:-1] = [var(indx) for indx in range(N)]
#    if isinstance(var, np.ndarray):
#        out = np.empty(N + 2, dtype=var.dtype)
#        out[1:-1] = var
#    return out


def get_dxvec(L, Nvol):
    if Nvol["a"]:
        dxa = Nvol["a"] * [L["a"]/Nvol["a"]]
    else:
        dxa = []
    if Nvol["s"]:
        dxs = Nvol["s"] * [L["s"]/Nvol["s"]]
    else:
        dxs = []
    dxc = Nvol["c"] * [L["c"]/Nvol["c"]]
    out = np.array(dxa + dxs + dxc)
    return out
