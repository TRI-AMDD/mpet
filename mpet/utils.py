import subprocess as subp

import numpy as np

import daetools.pyDAE as dae


def set_plot_defaults(mpl):
    """Set list of matplotlib rc parameters to make more readable plots."""
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
    """Calculate the linear mean along a vector or between two values."""
    if isinstance(a, np.ndarray):
        return 0.5*(a[1:] + a[:-1])
    else:
        return 0.5*(a + b)


def mean_harmonic(a, b=None):
    """Calculate the harmonic mean along a vector or between two values."""
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
    """Convert a constant to an array of length N."""
    out = np.array([val for indx in range(N)], dtype=object)
    return out


def get_var_vec(var, N, dt=False):
    """Convert a dae tools variable to a numpy array. Optionally return the time derivative of the
    variable.
    """
    if dt is True:
        out = np.array([var.dt(indx) for indx in range(N)])
    else:
        out = np.array([var(indx) for indx in range(N)])
    return out


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
    return out


def get_dxvec(L, Nvol):
    """Get a vector of cell widths spanning the full cell."""
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


def get_git_info(local_dir, shell=False):
    commit_hash = subp.check_output(['git', '-C', local_dir, 'rev-parse', '--short', 'HEAD'],
                                    stderr=subp.STDOUT, universal_newlines=True, shell=shell)
    commit_diff = subp.check_output(['git', '-C', local_dir, 'diff'],
                                    stderr=subp.STDOUT, universal_newlines=True, shell=shell)
    branch_name = subp.check_output(
        ['git', '-C', local_dir, 'rev-parse', '--abbrev-ref', 'HEAD'],
        stderr=subp.STDOUT, universal_newlines=True)
    return branch_name, commit_hash, commit_diff
