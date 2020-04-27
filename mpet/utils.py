import subprocess as subp

import numpy as np

import daetools.pyDAE as dae


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


def get_negative_sign_change_arrays(input_array):
    """This function takes an array of (+1, +1, +1, -1, -1, -1... +1, -1 ...) and splits it
       into a number of arrays in the y direction which are (0, 0, 0, 1, 1, 1... 0, 0) 
       whenever the array hits a sign change. It should have the number of cycles as rows.
       Thus it will be size (N*M) for each output, where N is the number of cycles 
       and M is the size of the array. In each ith row there are only 1's for the ith charging
       cycle.
       pos_segments is an array of charging segments, and neg_segments
       is true if we want arrays of discharge segments.
       """
    sign_mults = np.zeros((len(input_array) - 1)) #is +1 if no sign change, -1 if sign change@i+1
    for i in range(len(input_array)-1):
        sign_mults[i] = (input_array[i] * input_array[i+1] - 1) / (-2)
    #ends up 0 if no sign change, +1 if sign change
    #the odd sign changes indicate the beginning of the discharge cycle
    indices = np.array(np.nonzero(sign_mults)).flatten() #get indices of nonzero elements
    neg_indices_start = indices[::2] + 1
    neg_indices_end = indices[1::2] + 1
    pos_indices_start = indices[1::2] + 1
    pos_indices_start = np.delete(pos_indices_start, -1)
    pos_indices_start = np.insert(pos_indices_start, 0, 0)
    pos_indices_end = indices[::2] + 1
    #gets the beginning and ends of the charge and discharge cycles
    #pos_start, pos_end, neg_start, neg_end
    #if we did not finish the simulation, then we will 
    tot_cycle_number = len(pos_indices_start)
    #starts counting charge cycle at beginning and removes rest cycle element
    output_array_neg = np.zeros((tot_cycle_number, len(input_array)))
    output_array_pos = np.zeros((tot_cycle_number, len(input_array)))
    for j in range(tot_cycle_number):
    #for each segment
       #discharge segment
       ones_array = np.ones(neg_indices_end[j] - neg_indices_start[j])
       output_array_neg[j, neg_indices_start[j]:neg_indices_end[j]] = ones_array
       #fills with ones if it is the discharge segment number j
       #charge segment
       ones_array = np.ones(pos_indices_end[j] - pos_indices_start[j])
       output_array_pos[j, pos_indices_start[j]:pos_indices_end[j]] = ones_array
       #fills with ones if it is the charge segment number j
    return output_array_neg, output_array_pos


def get_density(material_type):
    """Gets active material density from input material type, in units of kg/m^3"""
    if material_type == "LiMn2O4": #cathode, so do LiMn2O4
        return 4.01e3
    elif material_type == "LiC6":
        return 2.26e3  #anode, so do graphite
    elif material_type == "NCA":
        return 4.45e3
    elif material_type == "LFP":
        return 3.6e3
        #https://cdn.intechopen.com/pdfs/18671/InTech-Lifepo4_cathode_material.pdf
#
#
#def get_mol_weight(material_type):
#    """Gets active material weight from input material type"""
#    if material_type == "LiMn2O4":
#        return 180.8
#    elif material_type == "LiC6":
#        return 12.01
#    elif material_type == "NCA":
#        return 183.54
#    elif material_type == "LFP":
#        return 157.757
