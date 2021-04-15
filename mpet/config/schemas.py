# definition of config file sections, parameters, types

import ast
from distutils.util import strtobool

from schema import Schema, Use, Optional, And, Or
import numpy as np


def parse_segments(key):
    """
    Parse the segments key of the configuration file and
    validate it
    """
    segments = ast.literal_eval(key)
    assert isinstance(segments, list), "segments must be a list"
    assert len(segments) > 0, "There must be at least one segment"
    for item in segments:
        assert (isinstance(item, tuple)) and (len(item) == 2), \
            "Each segment must be a tuple of (setpoint, time)"
    return segments


def check_allowed_values(value, allowed_values):
    """
    Check if value was chosen from a set of allowed values
    """
    assert value in allowed_values, f"{value} is invalid, options are: {allowed_values}"
    return True


def tobool(value):
    """
    Convert value to boolean
    """
    return bool(strtobool(value))


# system parameters, per section
system = {'Sim Params': {'profileType': lambda x:
                         check_allowed_values(x, ["CC", "CV", "CCsegments", "CVsegments"]),
                         'Crate': Use(float),
                         Optional('1C_current_density', default=None): Use(float),
                         Optional('tramp', default=0.): Use(float),
                         'Vmax': Use(float),
                         'Vmin': Use(float),
                         'Vset': Use(float),
                         Optional('capFrac', default=1.0): Use(float),
                         'segments': Use(parse_segments),
                         Optional('prevDir', default=''): str,
                         'tend': And(Use(float), lambda x: x > 0),
                         'tsteps': And(Use(int), lambda x: x > 0),
                         'relTol': And(Use(float), lambda x: x > 0),
                         'absTol': And(Use(float), lambda x: x > 0),
                         'T': Use(float),
                         'randomSeed': Use(tobool),
                         Optional('seed'): And(Use(int), lambda x: x >= 0),
                         Optional('dataReporter', default='mat'): str,
                         'Rser': Use(float),
                         'Nvol_c': And(Use(int), lambda x: x > 0),
                         'Nvol_s': And(Use(int), lambda x: x >= 0),
                         'Nvol_a': And(Use(int), lambda x: x >= 0),
                         'Npart_c': And(Use(int), lambda x: x >= 0),
                         'Npart_a': And(Use(int), lambda x: x >= 0)},
          'Electrodes': {'cathode': str,
                         'anode': str,
                         'k0_foil': Use(float),
                         'Rfilm_foil': Use(float)},
          'Particles': {'mean_c': Use(float),
                        'stddev_c': Use(float),
                        'mean_a': Use(float),
                        'stddev_a': Use(float),
                        'cs0_c': Use(float),
                        'cs0_a': Use(float),
                        Optional('specified_psd_c', default=False):
                            Or(Use(tobool), Use(lambda x: np.array(ast.literal_eval(x)))),
                        Optional('specified_psd_a', default=False):
                            Or(Use(tobool), Use(lambda x: np.array(ast.literal_eval(x))))},
          'Conductivity': {'simBulkCond_c': Use(tobool),
                           'simBulkCond_a': Use(tobool),
                           'sigma_s_c': Use(float),
                           'sigma_s_a': Use(float),
                           'simPartCond_c': Use(tobool),
                           'simPartCond_a': Use(tobool),
                           'G_mean_c': Use(float),
                           'G_stddev_c': Use(float),
                           'G_mean_a': Use(float),
                           'G_stddev_a': Use(float)},
          'Geometry': {'L_c': Use(float),
                       'L_a': Use(float),
                       'L_s': Use(float),
                       'P_L_c': Use(float),
                       'P_L_a': Use(float),
                       'poros_c': Use(float),
                       'poros_a': Use(float),
                       'poros_s': Use(float),
                       'BruggExp_c': Use(float),
                       'BruggExp_a': Use(float),
                       'BruggExp_s': Use(float)},
          'Electrolyte': {'c0': Use(float),
                          'zp': Use(int),
                          'zm': And(Use(int), lambda x: x < 0),
                          'nup': Use(int),
                          'num': Use(int),
                          'elyteModelType': str,
                          'SMset': str,
                          'n': Use(int),
                          'sp': Use(int),
                          'Dp': Use(float),
                          'Dm': Use(float)}}

electrode = {'Particles': {'type': str,
                           'discretization': Use(float),
                           'shape': str,
                           Optional('thickness'): Use(float)},
             'Material': {'muRfunc': str,
                          'logPad': Use(tobool),
                          'noise': Use(tobool),
                          'noise_prefac': Use(float),
                          'numnoise': Use(int),
                          Optional('Omega_a', default=None): Use(float),
                          Optional('Omega_b', default=None): Use(float),
                          Optional('Omega_c', default=None): Use(float),
                          'kappa': Use(float),
                          'B': Use(float),
                          Optional('EvdW', default=None): Use(float),
                          'rho_s': Use(float),
                          'D': Use(float),
                          'Dfunc': str,
                          Optional('E_D', default=0.): Use(float),
                          'dgammadc': Use(float),
                          'cwet': Use(float)},
             'Reactions': {'rxnType': str,
                           'k0': Use(float),
                           Optional('E_A', default=0.): Use(float),
                           'alpha': Use(float),
                           'lambda': Use(float),
                           'Rfilm': Use(float)}}


# convert the dictionaries to actual schemas
for d in [system, electrode]:
    for key, value in d.items():
        d[key] = Schema(value, ignore_extra_keys=False)
