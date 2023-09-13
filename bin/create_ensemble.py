#!/usr/bin/env python3

import configparser
import sys
import itertools
import os


def ensemble_definitions():
    # Values that need ensemble
    ensemble = [
        [("Geometry","L_c"), ["2","43"]],
        [("Geometry","L_a"), ["3","5","43"]],
    ]

    # helpers
    keys = [vals[0] for vals in ensemble]
    val = [vals[1] for vals in ensemble]
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
