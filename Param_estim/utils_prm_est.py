# utils_prm_est.py
import random
import os
from instructions import storage_folder


def mat_param_name(ensambles_tot_new):
    # create a name from the parameters
    nicename_mat = []
    for ensambles in ensambles_tot_new:
        # if anode is not present
        if len(ensambles) == 0:
            continue
        for val in ensambles:
            # convert val[1][0] in float, remove excess numbers after comma 
            # and convert in string
            float_val = float(val[1][0])
            float_val = "{:.4e}".format(float_val)
            value_name = str(float_val)
            nicename_mat.append(val[0][1] + "=" + str(value_name))
    nicename_mat = "_".join(nicename_mat)
    return nicename_mat

def make_folders(ensambles_tot_new, temperatures):
    # create store folder
    if not os.path.exists(storage_folder()):
        os.makedirs(storage_folder())

    store_fold_mat = mat_param_name(ensambles_tot_new)
    store_fold_mat = os.path.join(storage_folder(), store_fold_mat)
    # create store folder for that parameter set
    if not os.path.exists(store_fold_mat):
        os.makedirs(store_fold_mat)

    for temps in temperatures:
        store_folde_temp = os.path.join(store_fold_mat, "T=" + str(temps))
        # create store folder for that temperature
        if not os.path.exists(store_folde_temp):
            os.makedirs(store_folde_temp)
    return store_fold_mat

def generate_tuple_bounds(ensamble):
    parameter_tuples = []
    for ensables in ensamble:
        for param_tuple in ensables:
            if len(param_tuple[1]) == 2:
                lower_limit, upper_limit = param_tuple[1]
                parameter_tuples.append((float(lower_limit), float(upper_limit)))
            else:
                raise Exception('Parameter range has to be of length 2')

    return tuple(parameter_tuples)


def generate_initial_guess(tuple_bounds, rand = True, seed = 1):
    random.seed(seed)
    initial_guess = []
    for bounds in tuple_bounds:
        if rand:
            initial_guess.append(random.uniform(bounds[0], bounds[1]))
        else:
            initial_guess.append((bounds[0] + bounds[1])/2)
    return initial_guess