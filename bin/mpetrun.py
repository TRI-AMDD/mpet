#!/usr/bin/env python3

import sys
import argparse
from argparse import RawTextHelpFormatter

from mpet.version import __version__
import mpet.main as main

desc = """MPET - Multiphase Porous Electrode Theory
This software is designed to run simulations of batteries with porous electrodes
using porous electrode theory, which is a volume-averaged, multiscale approach
to capture the coupled behavior of electrolyte and active material within
electrodes.

If you use this software in academic work, please cite the relevant references
detailing its development as presented in the LICENSE file.

See also: https://bitbucket.org/bazantgroup/mpet"""

parser = argparse.ArgumentParser(description=desc, formatter_class=RawTextHelpFormatter)
parser.add_argument('file', help='MPET system configuration file')
parser.add_argument('-v','--version', action='version',
                    version='%(prog)s '+__version__)
args = parser.parse_args()

try:
    main.main(sys.argv[1])
except IndexError:
    print("ERROR: No parameter file specified. Aborting")
    raise
