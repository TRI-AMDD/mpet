# MPET -- Multiphase Porous Electrode Theory
This software is design to run simulations of batteries with porous electrodes using porous
electrode theory, which is a volume-averaged, multiscale approach to capturing the coupled behavior
of electrolyte and active material within electrodes. As a result, with physical paramter intputs
and run protocols (specified current or voltage profiles), it makes predictions about the internal
dynamics within a battery (electrolyte concentration and potential, solid phase concentrations,
reaction rates, etc.) and also macroscopic, easily measurable electrochemical quantities such as
total current and voltage. In this way, it is similar the
[dualfoil](http://www.cchem.berkeley.edu/jsngrp/fortran.html) code released by Newman and coworkers
from Berkeley. This software has most of the functionality contained in dualfoil (it is currently
missing temperature dependence). However, beyond the standard porous electrode theory
simulations, this software can also simulate electrodes in which the active materials phase
separate using non-equilibrium thermodynamics within a phase field modeling framework. Such
behavior is common in widely used electrode materials, including for example graphite and LiFePO4.

## Installation

1.  Install `python3.4`, `numpy`, `scipy`, `matplotlib`, `pyqt4`
    - Linux: Try using your package manager
    - Windows: Use Anaconda
        - Remove any other Python installation
        - Get and install [Anaconda](https://www.continuum.io/downloads) (32 bit, Python 3.5)
        - Open `cmd.exe`
        - `$ conda install python=3.4 anaconda` to convert to a Python 3.4 installation
        - Say "y" when prompted. This may take several minutes.
2.  Install [DAE Tools](https://sourceforge.net/projects/daetools/files/1.4.0)
    - Get the version corresponding to your operating system (and py34 if on Windows)
3.  Download a copy of this repository to some place on your system path (for example, put this
    directory within a working diretory in which you want to run simulations).

If you want to use DAE Tools with a different version of Python3, you can compile it from source as
described [here](http://daetools.com/docs/getting_daetools.html).

## Simulation

1.  Copy `mpetrun` and `mpetplot` from the repository to your working directory (the directory from
    which you'll run simulations and in which output will be stored).
2.  Copy the overall system parameters file,
    `configDefaults/params_system.cfg`, to your working directory .
3.  Copy at least one material parameters file from `configDefaults`
    (e.g. `configDefaults/params_electrodes.cfg`) to the working directory
4.  Edit `params_system.cfg` to suit the simulation you're trying to run. Be
    sure to reference a material parameters file for the cathode and
    optionally one (the same or separate file) for the anode.
5.  Edit the material parameters file(s) serving as the electrode
    materials.
6.  Run `mpetrun`, passing `params_system.cfg` as an argument:
    `mpetrun params_system.cfg` on Windows or `./mpetrun params_system.cfg` on Linux/Mac

The software will save the simulation output data in a folder called `sim_output` and will also
keep a time-stamped copy in a folder called `history`. Each output directory should contain
- the output data (`.mat` file)
- copies of the input parameters files defining the simulation
- a copy of the daetools config parameters (e.g. solver tolerances)
- information about the script used to run the simulation
- information about the simulation (e.g. run time)
- processed, dimensional and nondimensional parameters as
  Python-pickled dictionary objects

Note, if you are seeing errors about reaching the maximum number of steps with suggestions about
scaling your tolerances, try increasing the IDAS:MaxNumSteps parameter to 100000. This can be found
in the `daetools.cfg` file. In Linux, this file is found within `/etc/daetools`, and in Windows, it
is within `C:\daetools`.

## Analysis

1.  Analyze output with `mpetplot` (pass output data directory, then
    plot-type as arguments)
    - e.g., voltage plot: `mpetplot sim_output v` or `./mpetplot sim_output v`
    - other options (`full`, `c`, `a` indicate full cell, cathode, and anode):
      - `v` -- voltage vs filling fraction
      - `curr` -- current vs time
      - `elytec` -- electrolyte concentration (movie)
      - `elytep` -- electrolyte potential (movie)
      - `csld_{c,a}` -- solid concentrations (all, movie, used with `solidType_{c,a}` not homog)
      - `phisld_{c,a}` -- solid potential (all, movie, used with `simSurfCond_{c,a}` = true)
      - `cbar_{full,c,a}` -- average solid concentrations (movie)
      - `bulkp_{c,a}` -- macroscopic electrode solid phase potential (movie)
      - `soc_{c,a}` -- electrode state of charge
2.  Alternatively, convert the output to plain text (csv) format using
    `mpetplot sim_output text` (or replace `sim_output` with
    any subfolder in the `history` folder). Then analyze using whatever
    tools you prefer.

If you want to save output to a movie (or figure), add `save` as an extra
argument to `mpetplot`: `mpetplot sim_output cbar save`.

Movie output requires that you have `ffmpeg` or `mencoder` (part of
`MPlayer`) installed.

## Testing

When adding new features or making changes to the code, it's helpful
to run a suite of tests to make sure various things are behaving as
expected. This should not be necessary for a user who is not changing
the code at all, although it could still be nice to verify that the
outputs you are seeing match those the developer(s) expect for a few
specific cases. To run the tests do the following:

1. Download `ref_outputs.zip` from
   [here](http://mit.edu/smithrb/www/ref_outputs.zip) and unzip within
   the tests subdirectory. This should give make a number of
   directories of the form `tests/ref_outputs/test###`.
2. Copy `mpettest` to the working directory and run it. This will
   run a number of simulations and compare their output to those
   outputs from the downloaded reference outputs along with a few
   comparisons of simulations to their corresponding analytical
   results.
