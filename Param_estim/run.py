import os
from scipy.optimize import differential_evolution, minimize
from utils import generate_initial_guess, generate_tuple_bounds, scale
from objective import obj
from instructions import *


ensamble_tot, operating_conditions, config_files, optimization_method = sim_instructions()
if not os.path.exists(save_folder()):
    os.makedirs(save_folder())

with open(os.path.join(save_folder(),"log_book.txt"), 'w') as f:
    for ensambles in ensamble_tot:
        for val in ensambles:
            f.write(val[0][1] + "\t")
    f.write("mse\n")

scaled_ensamble, scaling = scale(ensamble_tot)
scaled_bounds = generate_tuple_bounds(scaled_ensamble)
scaled_ig = generate_initial_guess(scaled_bounds, rand=True, seed=0)

if optimization_method == 'DE':
    # Set options for DE
    options = {
        'strategy': 'best1bin',  # Strategy for DE
        'popsize': 10,           # Population size
        'tol': 1e-4,             # Tolerance for convergence
        'maxiter': 500,         # Maximum number of iterations
        'disp': False,           # Display convergence messages
        'seed': 0               # Random seed for reproducibility
    }

    # Run Differential Evolution optimization
    result = differential_evolution(obj, scaled_bounds, seed=options['seed'], strategy=options['strategy'],
                                    popsize=options['popsize'], tol=options['tol'], maxiter=options['maxiter'],
                                    disp=options['disp'])

elif optimization_method == 'GD':
    # Run Gradient descent

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
    }

    solution = minimize(obj, scaled_ig, method='L-BFGS-B', bounds=scaled_bounds, options=options)
else:
    print('Optimization method not recognized')
    print('Available methods are: DE, GD')

with open(os.path.join(save_folder(),"log_book.txt"), 'a') as f:
    f.write('final parameters' + "\n")
    for val in result.x:
        f.write(str(val) + "\t")
    f.write("\n")
    f.write("end\n")
