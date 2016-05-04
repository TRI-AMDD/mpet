# MPET

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

If you want to use DAE Tools with a different version of Python3, you can compile it from source as
described [here](http://daetools.com/docs/getting_daetools.html).

## Simulation

1. Copy the overall system parameters file,
   `configDefaults/params_system.cfg`, to main directory
2. Copy at least one material parameters file from `configDefaults`
   (e.g. `configDefaults/params_electrodes.cfg`) to the main directory
3. Edit `params_system.cfg` to suit the simulation you're trying to run. Be
   sure to reference a material parameters file for the cathode and
   optionally one (the same or separate file) for the anode.
4. Edit the material parameters file(s) serving as the electrode
   materials
5. Run script, passing `params_system.py` as an argument (try with/without the 3, or with 3.4):
   `python[3] mpet.py params_system.py`

The software will save the simulation output data in a folder called `sim_output` and will also
keep a time-stamped copy in a folder called `history`. Each output directory should contain
- the output data (`.mat` file)
- copies of the input parameters files defining the simulation
- a copy of the daetools config parameters (e.g. solver tolerances)
- information about the script used to run the simulation
- information about the simulation (e.g. run time)
- processed, dimensional and nondimensional parameters as
  Python-pickled dictionary objects

Note, if you are seeing errors about reaching the maximum number of
steps with suggestions about scaling your tolerances, try increasing
the IDAS:MaxNumSteps parameter to 100000. This can be found in the
`daetools.cfg` file. In Linux, this file is found within
`/etc/daetools`, and in Windows, it is within `C:\daetools`.

## Analysis

1.  Analyze output with `plot_data.py` (pass output data directory, then
    plot-type as arguments)
    - e.g., voltage plot: `python plot_data.py sim_output v`
    - other options (`full`, `c`, `a` indicate full cell, cathode, and anode):
      - `v` -- voltage vs filling fraction
      - `curr` -- current vs time
      - `elytec` -- electrolyte concentration (movie)
      - `elytep` -- electrolyte potential (movie)
      - `csld_{c,a}` -- solid concentrations (all, movie, used with `solidType_{c,a}` not homog)
      - `phisld_{c,a}` -- solid potential (all, movie, used with `simSurfCond_{c,a}` = true)
      - `cbar_{full,c,a}` -- average solid concentrations (movie)
      - `bulkp_{c,a}` -- macroscopic electrode solid phise potential (movie)
      - `soc_{c,a}` -- electrode state of charge
2.  Alternatively, convert the output to plain text (csv) format using
    `python outmat2txt.py sim_output` (or replace `sim_output` with
    any subfolder in the `history` folder). Then analyze using whatever
    tools you prefer.

If you want to save output to a movie (or figure), add `save` as an extra
argument to `plot_data.py`: `python plot_data.py sim_output cbar save`.

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
2. Within the tests directory, run `python test_suite.py`. This will
   run a number of simulations and compare their output to those
   outputs from the downloaded reference outputs along with a few
   comparisons of simulations to their corresponding analytical
   results.
