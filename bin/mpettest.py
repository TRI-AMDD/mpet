#!/usr/bin/env python3

import errno
import os
import os.path as osp
import sys
import time

import tests.test_suite as tst


def silentremove(filename):
    try:
        os.remove(filename)
    except OSError:
        pass

 
if len(sys.argv) > 1:
    compareDir = osp.join(os.getcwd(), sys.argv[1])
else:
    compareDir = None
timeStart = time.time()
tst.main(compareDir)
timeEnd = time.time()
tTot = timeEnd - timeStart
print("Total test time:", tTot, "s")

#remove files if they exist
 
silentremove("LA4_8rep.MWF") 
silentremove("Short_PreDiag_000173.000") 
silentremove("test_mwf_LA4.000")

