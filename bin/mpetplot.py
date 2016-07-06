#!/usr/bin/env python3

import os
import sys

import matplotlib as mpl
# DAE Tools uses Qt4Agg, so we might as well use the same.
# Any installed interactive backend other than Qt5Agg should work.
mpl.use("Qt4Agg")
import matplotlib.pyplot as plt

import mpet.outmat2txt as outmat2txt
import mpet.plot_data as plot_data

# Get input file from script parameters
if len(sys.argv) < 2:
    raise Exception("Need input data directory name")
indir = sys.argv[1]
if not os.path.exists(os.path.join(os.getcwd(), indir)):
    raise Exception("Input file doesn't exist")
# Optionally just convert output to text
if len(sys.argv) == 3 and sys.argv[2] == "text":
    outmat2txt.main(indir)
    sys.exit()
# Get plot type from script parameters
plots = []
if len(sys.argv) > 2:
    plots.append(sys.argv[2])
else:
    plots.append("v")
# Save the plot instead of showing on screen?
# Get from script parameters
save_flag = False
print_flag = True
data_only = False
save_only = False
if len(sys.argv) > 3:
    if sys.argv[3] in ["save", "saveonly"]:
        save_flag = True
        if sys.argv[3] == "saveonly":
            save_only = True
            print_flag = False
    else:
        for i in range(3, len(sys.argv)):
            plots.append(sys.argv[i])
out = []
for plot_type in plots:
    out.append(plot_data.show_data(
        indir, plot_type, print_flag, save_flag, data_only))
if not save_only:
    plt.show()
