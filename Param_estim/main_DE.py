import os
from scipy.optimize import differential_evolution
from utils_prm_est import generate_initial_guess, generate_tuple_bounds
from objective import obj
from instructions import *
from util_run import scale

# Generate initial guesses for all parameters
ensamble_tot, operating_conditions, config_files = sim_instructions()

# Take names of parameters from ensamble_tot
# Open log_book.txt in storage folder
# Write the name of the parameters
store_folder = storage_folder()
    # create store folder
if not os.path.exists(store_folder):
    os.makedirs(store_folder)
with open(os.path.join(storage_folder(),"log_book.txt"), 'w') as f:
    for ensambles in ensamble_tot:
        for val in ensambles:
            f.write(val[0][1] + "\t")
    f.write("mse\n")

scaled_ensamble, scaling = scale(ensamble_tot)
scaled_bounds = generate_tuple_bounds(scaled_ensamble)
scaled_ig = generate_initial_guess(scaled_bounds, rand=True, seed=0)


# Set options for DE
options = {
    'strategy': 'best1bin',  # Strategy for DE
    'popsize': 12,           # Population size
    'tol': 1e-6,             # Tolerance for convergence
    'maxiter': 1000,         # Maximum number of iterations
    'disp': False,           # Display convergence messages
    'seed': 0                # Random seed for reproducibility
}

# Run Differential Evolution optimization
result = differential_evolution(obj, scaled_bounds, seed=options['seed'], strategy=options['strategy'],
                                 popsize=options['popsize'], tol=options['tol'], maxiter=options['maxiter'],
                                 disp=options['disp'])



with open(os.path.join(storage_folder(),"log_book.txt"), 'a') as f:
    f.write('final parameters' + "\n")
    for val in result.x:
        f.write(str(val) + "\t")
    f.write("\n")
    f.write("end\n")
