import subprocess as subp

import os
import numpy as np
import h5py
import re
import scipy.io as sio

import daetools.pyDAE as dae

from sympy.parsing.sympy_parser import parse_expr


def mean_linear(a):
    """Calculate the linear mean along a vector."""
    return 0.5*(a[1:] + a[:-1])


def weighted_linear_mean(a, wt):
    return (wt[1:]*a[1:] + wt[:-1]*a[:-1])/(wt[1:]+wt[:-1])


def mean_harmonic(a):
    """Calculate the harmonic mean along a vector."""
    return (2 * a[1:] * a[:-1]) / (a[1:] + a[:-1] + 1e-20)


def weighted_harmonic_mean(a, wt):
    return((wt[1:]+wt[:-1])/(wt[1:]/a[1:]+wt[:-1]/a[:-1]))


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
            try:
                varout[sectn] = np.zeros(Nvol[sectn])
            except KeyError:
                varout[sectn] = np.zeros(0)
    out = np.hstack((varout["a"], varout["s"], varout["c"]))
    return out


def get_dxvec(L, Nvol):
    """Get a vector of cell widths spanning the full cell."""
    if "a" in Nvol:
        dxa = Nvol["a"] * [L["a"]/Nvol["a"]]
    else:
        dxa = []
    if "s" in Nvol:
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


def open_data_file(dataFile):
    """Load hdf5/mat file output.
    Always defaults to .mat file, else opens .hdf5 file.
    Takes in dataFile (path of file without .mat or .hdf5), returns data (output of array)"""
    data = []
    if os.path.isfile(dataFile + ".mat"):
        data = sio.loadmat(dataFile + ".mat")
    elif os.path.isfile(dataFile + ".hdf5"):
        data = h5py.File(dataFile + ".hdf5", 'r')
    else:
        raise Exception("Data output file not found for either mat or hdf5 in " + dataFile)
    return data


def get_dict_key(data, string, squeeze=True, final=False):
    """Gets the values in a 1D array, which is formatted slightly differently
    depending on whether it is a h5py file or a mat file
    Takes in data array and the string whose value we want to get. Final is a
    boolean that determines whether or not it only returns the final value of the array.
    If final is true, then it only returns the last value, otherwise it returns the entire array.
    Final overwrites squeeze--if final is true, then the array will always be squeezed.
    Squeeze squeezes into 1D array if is true, otherwise false"""
    # do not call both squeeze false and final true!!!
    if final:  # only returns last value
        return data[string][...,-1].item()
    elif squeeze:
        return np.squeeze(data[string][...])
    else:  # returns entire array
        return data[string][...]


def get_crate(crate, Cratecurr):
    """Returns crate from Crate input in params.sys.
    if it is in Crate, returns original number. otherwise converts
    from A to Crate and returns that number.
    if it is in waveform, returns waveform
    Cratecurr is theoretical 1C current density"""
    out = 0
    if "t" in str(crate):
        # if waveform, we output original formula
        if "A" not in str(crate)[-1]:
            out = parse_expr(crate)
        else:  # convert A to Crate
            amp_value = re.sub(r'[A]+', '', crate, re.I)
            amp_value = parse_expr(amp_value)
            out = amp_value / Cratecurr
    else:
        # if not waveform, output float
        if "A" not in str(crate)[-1]:
            out = float(crate)
        else:  # convert A to Crate
            amp_value = float(re.sub(r'[A]+', '', crate, re.I))
            out = amp_value / Cratecurr
    return out


def get_negative_sign_change_arrays(input_array):
    """This function takes an array of (+1, +1, +1, -1, -1, -1... +1, -1 ...) and splits it
       into a number of arrays in the y direction which are (0, 0, 0, 1, 1, 1... 0, 0)
       whenever the array hits a sign change. It should have the number of cycles as rows.
       Thus it will be size (N*M) for each output, where N is the number of cycles
       and M is the size of the array. In each ith row there are only 1's for the ith charging
       cycle.
       Returns beginning and end discharge and charge segments in order
       """
    sign_mults = np.zeros((len(input_array) - 1))  # is +1 if no sign change, -1 if sign change@i+1
    for i in range(len(input_array)-1):
        # ends up 0 if no sign change, +1 if sign change
        sign_mults[i] = (input_array[i] * input_array[i+1] - 1) / (-2)
    # if we have no outputs with sign change, then end
    if np.all(sign_mults == 0):
        print("ERROR: Did not complete a single cycle, cannot plot cycling plots")
        raise
    # the odd sign changes indicate the beginning of the discharge cycle
    indices = np.array(np.nonzero(sign_mults)).flatten()  # get indices of nonzero elements
    neg_indices_start = indices[::2] + 1
    neg_indices_end = indices[1::2] + 1
    pos_indices_start = indices[1::2] + 1
    pos_indices_start = np.delete(pos_indices_start, -1)
    pos_indices_start = np.insert(pos_indices_start, 0, 0)
    pos_indices_end = indices[::2] + 1
    return neg_indices_start, neg_indices_end, pos_indices_start, pos_indices_end


def get_density(material_type):
    """Gets active material density from input material type, in units of kg/m^3"""
    if material_type == "LMO":  # cathode, so do LiMn2O4
        return 4.01e3
    elif material_type == "LiC6":
        return 2.26e3  # anode, so do graphite
    elif material_type == "NCA":
        return 4.45e3
    elif material_type == "LFP":
        return 3.6e3
        # https://cdn.intechopen.com/pdfs/18671/InTech-Lifepo4_cathode_material.pdf
