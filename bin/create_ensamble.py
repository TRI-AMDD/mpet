#!/usr/bin/env python3

import configparser
import sys
import itertools
from copy import deepcopy
import os

#Values that need ensamble
ensamble = [
  [("Geometry","L_c")  , ["2","43"]],
  [("Geometry","L_a")  , ["3","5","43"]],
]

#helpers
keys = [ vals[0] for vals in ensamble ]
val  = [ vals[1] for vals in ensamble ]


if __name__ == '__main__':
  # Read in file
  if len(sys.argv) < 2:
    print("need the config file [python create_enamble.py <baseconfig>]")
    exit(1)
  cfg = configparser.ConfigParser()
  cfg.optionxform = str
  cfg.read(sys.argv[1])

  #Create all variations
  combinations = list(itertools.product(*val))

  for combination in combinations:
    params = dict(zip(keys,combination))
    print(params)

    new_cfg = cfg
    
    nicename = []
    for key, val in params.items():
      new_cfg[key[0]][key[1]] = val
      nicename.append(key[1] + "=" + val)


    # Write config 
    cfg_dir = os.path.dirname(sys.argv[1])
    with open(cfg_dir + "/" + "-".join(nicename) + ".cfg", "w") as f:
      new_cfg.write(f)
