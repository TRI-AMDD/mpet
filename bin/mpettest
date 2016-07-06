#!/usr/bin/env python3

import os
import os.path as osp
import sys
import time

# DAE Tools uses Qt4Agg, so we might as well use the same.
# Any installed interactive backend other that Qt5Agg should work.
import matplotlib
matplotlib.use("Qt4Agg")

import mpet.tests.test_suite as tst

if len(sys.argv) > 1:
    compareDir = osp.join(os.getcwd(), sys.argv[1])
else:
    compareDir = None
timeStart = time.time()
tst.main(compareDir)
timeEnd = time.time()
tTot = timeEnd - timeStart
print("Total test time:", tTot, "s")
