#!/usr/bin/env python3

import sys

import mpet.main as main

try:
    main.main(sys.argv[1])
except IndexError:
    print("ERROR: No parameter file specified. Aborting")
    raise
