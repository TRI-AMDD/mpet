# main.py

from scipy.optimize import minimize
from utils_prm_est import generate_initial_guess, generate_tuple_bounds
from objective import obj
from instructions import sim_instructions, storage_folder
import os
from util_run import scale

# create txt file containing paramters and mse in strage folders

# Generate initial guesses for all parameters
# Possible to select random guess within the bounds
# or to use the mean of the bounds
ensamble_tot, operating_conditions, config_files = sim_instructions()

# take name of paramters from ensamble_tot
# open log_book.txt in storage folder
# write the name of the paramters

with open(os.path.join(storage_folder(), "log_book.txt"), 'w') as f:
    for ensambles in ensamble_tot:
        for val in ensambles:
            f.write(val[0][1] + "\t")
    f.write("mse\n")

scaled_ensamble, scaling = scale(ensamble_tot)
scaled_bounds = generate_tuple_bounds(scaled_ensamble)
scaled_ig = generate_initial_guess(scaled_bounds, rand=True, seed=0)

options = {
    'maxcor': 50,       # Maximum number of variable metric corrections (limited-memory BFGS)
    # 'ftol': 2.2e-9, 
    'ftol': 1e-8,     # Tolerance for the relative change in the objective function
    'gtol': 1e-6,       # Tolerance for the relative change in the gradient of the objective function
    'eps': 1e-5,        # Small positive number to avoid division by zero in gradient computations
    'maxfun': 1500,    # Maximum number of function evaluations allowed
    'maxiter': 1500,   # Maximum number of iterations allowed
    'iprint': -1,       # Controls the verbosity of the output (0: no output, 1: minimal, 2: detailed)
    'disp': False,      # If True, displays convergence messages during optimization
    'maxls': 20,        # Maximum number of line search steps in each iteration
    'ftolmin': 1e-12,
    # 'ftolmin': 1e-12,   # Minimum tolerance for convergence in the objective function value
    'maxstep': None     # Maximum step length in a line search (None: unconstrained)
}

solution = minimize(obj, scaled_ig, method='L-BFGS-B', bounds=scaled_bounds, options=options)

with open(os.path.join(storage_folder(), "log_book.txt"), 'a') as f:
    f.write('final parameters' + "\n")
    for val in solution.x:
        f.write(str(val) + "\t")
    f.write("\n")
    f.write("end\n")