# MPET -- Multiphase Porous Electrode Theory

This software is designed to run simulations of batteries with porous electrodes using porous
electrode theory, which is a volume-averaged, multiscale approach to capture the coupled behavior
of electrolyte and active material within electrodes. As a result, with physical parameter inputs
and run protocols (specified current or voltage profiles), it makes predictions about the internal
dynamics within a battery (electrolyte concentration and potential, solid phase concentrations,
reaction rates, etc.) and also macroscopic, easily measurable electrochemical quantities such as
total current and voltage. In this way, it is similar to the
[`dualfoil`](http://www.cchem.berkeley.edu/jsngrp/fortran.html) code released by Newman and
coworkers from Berkeley. This software has much of the functionality contained in `dualfoil` (it is
currently missing, e.g., temperature dependence). However, beyond the standard porous electrode
theory simulations, this software can also simulate electrodes in which the active materials phase
separate using non-equilibrium thermodynamics within a phase field modeling framework. Such
behavior is common in widely used electrode materials, including graphite and LiFePO4.

## Installation

1.  Install `python3.4`, `numpy`, `scipy`, `matplotlib`, `pyqt4` (Note: `pyqt4` should be safe to
    neglect if you don't plan to run any DAE Tools tutorials.)
    - Linux: Try using your package manager or using Anaconda as in the Windows instructions.
    - Windows: Use Anaconda
        - Remove any other Python installation
        - Get and install [Anaconda](https://www.continuum.io/downloads) (32 bit, Python 3.x).
          Install for all users.
        - Open `cmd.exe`
        - `$ conda install python=3.4 anaconda` to convert to a Python 3.4 installation
        - Say "y" when prompted. This may take several minutes.
2.  Install the most recent version of [DAE Tools](https://sourceforge.net/projects/daetools/files)
    - Get the version corresponding to your operating system (and py34 if on Windows)
3.  Download a copy of this repository to some place on your system path (for example, put this
    directory within a working diretory in which you want to run simulations).

If you want to use DAE Tools with a different version of Python 3, you can compile it from source
as described [here](http://daetools.com/docs/getting_daetools.html).
Also note that DAE Tools can be installed within a python virtual environment, so feel free to take
that approach instead.

## Simulation

1.  Enter the root repository directory. This will serve as your working directory. Simulations
    will be run from here, and outputs will be stored here.
2.  Copy `mpetrun.py` and `mpetplot.py` from the `bin` directory to your working directory.
3.  Copy the overall system parameters file,
    `configs/params_system.cfg`, to your working directory.
3.  Copy at least one material parameters file from `configs`
    (e.g. `configs/params_electrodes.cfg`) to the working directory.
4.  Edit `params_system.cfg` to suit the simulation you're trying to run. Be
    sure to reference a material parameters file for the cathode and
    optionally one (the same or separate file) for the anode.
5.  Edit the material parameters file(s) serving as the electrode
    materials.
6.  Run `mpetrun.py`, passing `params_system.cfg` as an argument:
    `python mpetrun.py params_system.cfg`
    (or optionally `./mpetrun params_system.cfg` on Linux)

The software will save the simulation output in a time-stamped subdirectory within a directory
called `history`. The data contents of the most recent output will also be copied to a directory
called `sim_output`. Each output directory should contain
- the output data (`.mat` file)
- copies of the input parameters files defining the simulation
- a copy of the daetools config parameters (e.g. solver tolerances)
- information about the script used to run the simulation
- information about the simulation (e.g. run time)
- processed, dimensional and nondimensional parameters as
  Python-pickled dictionary objects

## Analysis

1.  Analyze output with `mpetplot.py` (pass output data directory, then
    plot-type as arguments)
    - e.g., voltage plot: `python mpetplot.py sim_output v`
    - other options (`full`, `c`, `a` indicate full cell, cathode, and anode):
      - `v` or `vt` -- voltage vs filling fraction or vs time
      - `curr` -- current vs time
      - `elytec{f}` -- electrolyte concentration (movie) or final snapshot with `f`
      - `elytep{f}` -- electrolyte potential (movie) or final snapshot with `f`
      - `elytei{f}` -- electrolyte current density (movie) or final snapshot with `f`
      - `surf_{c,a}` -- solid surface concentrations
      - `soc_{c,a}` -- overall utilization / state of charge of electrode
      - `csld_{c,a}` -- solid concentrations of particles in electrode (movie; used with `solidType_{c,a}` not homog)
      - `cbarLine_{c,a}` -- average concentration in each particle of electrode
      - `cbar_{full,c,a}` -- average solid concentrations as changing colors (movie)
      - `bulkp_{c,a}` -- macroscopic electrode solid phase potential (movie)
2.  Alternatively, convert the output to plain text (csv) format using
    `python mpetplot.py sim_output text` (or replace `sim_output` with
    any subfolder in the `history` folder). Then analyze using whatever
    tools you prefer.

If you want to save output to a movie (or figure), add `save` as an extra
argument to `mpetplot.py`: `python mpetplot.py sim_output cbar save`.

Movie output requires that you have `ffmpeg` or `mencoder` (part of
`MPlayer`) installed.

## Troubleshooting
If you are seeing errors about reaching the maximum number of steps with suggestions about
scaling your tolerances, try increasing the IDAS:MaxNumSteps parameter to 100000. This can be
found in the `daetools.cfg` file. This is found within the DAE Tools install directory.

Otherwise, please use the bitbucket website to file issues and/or bug reports with the software.

## Testing

When adding new features or making changes to the code, it's helpful
to run a suite of tests to make sure various things are behaving as
expected. This should not be necessary for users who are not changing
the code at all, although it could still be nice to verify that the
outputs users are seeing match those the developers expect for a few
specific cases. To run the tests, do the following:

1. Download `ref_outputs.zip` from
   [here](http://mit.edu/smithrb/www/ref_outputs.zip) and unzip within
   the tests subdirectory. This should give a number of
   directories of the form `tests/ref_outputs/test###`.
2. Copy `mpettest.py` to the working directory and run it. This will
   run a number of simulations and compare their outputs to those
   outputs from the downloaded reference outputs along with a few
   comparisons of simulations to their corresponding analytical
   results.

Note that tests may "fail" even when things are okay, resulting from small numerical differences.
If tests fail, it is helpful to look at the comparison plots generated by default within
`mpet/tests/test_outputs/[time-stamped-directory]/plots` to see if the differences seem
significant.
