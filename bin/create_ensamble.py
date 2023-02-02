#!/usr/bin/env python3

import configparser
import sys
import itertools
import os

# Values that need ensamble
# ensamble = [
#     [("Particles","mean_c"), ["100e-9","300e-9"]],
#     [("Particles","stddev_c"), ["10e-9","100e-9"]],
#     [("Sim Params","relTol"), ["1e-5","1e-7"]],
#     [("Sim Params","absTol"), ["1e-5","1e-7"]],
#     [("Electrodes", "cathode"), ['D=1e-18-k0=2.cfg','D=1e-18-k0=2.5.cfg','D=1e-18-k0=3.cfg',
#                                  'D=5e-18-k0=2.cfg','D=5e-18-k0=2.5.cfg','D=5e-18-k0=3.cfg',
#                                  'D=50e-18-k0=2.cfg','D=50e-18-k0=2.5.cfg','D=50e-18-k0=3.cfg',
#                                  'D=500e-18-k0=2.cfg','D=50e-18-k0=2.5.cfg','D=500e-18-k0=3.cfg']]
# ]

# ensamble = [
#     [("Material","D"), ["1e-18","5e-18","50e-18","500e-18"]],
#     [("Reactions","k0"), ["2","2.5","3"]],
# ]


def ensemble_definitions():
    # Values that need ensamble
    ensamble = [[("Particles","mean_c"), ["100e-9","300e-9"]],
                [("Particles","stddev_c"), ["10e-9","100e-9"]],
                [("Sim Params","relTol"), ["1e-5","1e-7"]],
                [("Sim Params","absTol"), ["1e-5","1e-7"]],
                [("Electrodes", "cathode"), ['D=1e-18-k0=2.cfg','D=1e-18-k0=2.5.cfg','D=1e-18-k0=3.cfg',
                                             'D=5e-18-k0=2.cfg','D=5e-18-k0=2.5.cfg','D=5e-18-k0=3.cfg',
                                             'D=50e-18-k0=2.cfg','D=50e-18-k0=2.5.cfg','D=50e-18-k0=3.cfg',
                                             'D=500e-18-k0=2.cfg','D=50e-18-k0=2.5.cfg','D=500e-18-k0=3.cfg']]
                ]

    # helpers
    keys = [vals[0] for vals in ensamble]
    val = [vals[1] for vals in ensamble]
    return keys, val


def create_ensemble(cff, keys=None, val=None):
    with open('ensemble_parallel_configs.txt', "w") as ff:
        if keys is None and val is None:
            keys, val = ensemble_definitions()
        cfg = configparser.ConfigParser()
        cfg.optionxform = str
        cfg.read(cff)
        # Create all variations
        combinations = list(itertools.product(*val))
        for combination in combinations:
            params = dict(zip(keys,combination))
            new_cfg = cfg
            nicename = []
            for key, val in params.items():
                new_cfg[key[0]][key[1]] = val
                nicename.append(key[1] + "=" + val)

            # Write config
            cfg_dir = os.path.dirname(cff)
            with open(cfg_dir + "/" + "-".join(nicename) + ".cfg", "w") as f:
                new_cfg.write(f)
                ff.write(str(cfg_dir + "/" + "-".join(nicename) + ".cfg\n"))
    return


if __name__ == '__main__':
    # Read in file
    if len(sys.argv) < 2:
        print("need the config file [python create_enamble.py <baseconfig>]")
        exit(1)
    create_ensemble(sys.argv[1])