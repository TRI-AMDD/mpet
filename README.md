Development: [![Coverage Status](https://coveralls.io/repos/github/TRI-AMDD/mpet-dev/badge.svg?branch=development)](https://coveralls.io/github/TRI-AMDD/mpet-dev?branch=development)
Master: [![Coverage Status](https://coveralls.io/repos/github/TRI-AMDD/mpet-dev/badge.svg?branch=master)](https://coveralls.io/github/TRI-AMDD/mpet-dev?branch=master)
# MPET -- Multiphase Porous Electrode Theory

This software is designed to run simulations of batteries with porous electrodes using porous electrode theory, which is a volume-averaged, multiscale approach to capture the coupled behavior of electrolyte and active material within electrodes. As a result, with physical parameter inputs and run protocols (specified current or voltage profiles), it makes predictions about the internal dynamics within a battery (electrolyte concentration and potential, solid phase concentrations, reaction rates, etc.) and also macroscopic, easily measurable electrochemical quantities such as total current and voltage. In this way, it is similar to the [`dualfoil`](http://www.cchem.berkeley.edu/jsngrp/fortran.html) code released by Newman and coworkers from Berkeley. This software has much of the functionality contained in `dualfoil` (it is currently missing, e.g., temperature dependence). However, beyond the standard porous electrode theory simulations, this software can also simulate electrodes in which the active materials phase separate using non-equilibrium thermodynamics within a phase field modeling framework. Such behavior is common in widely used electrode materials, including graphite and LiFePO4.

If you use this software in academic work, please cite the relevant references detailing its development as presented in the `LICENSE` file. For more details on the theory implemeneted in the code, see:

Smith, R. B., and Bazant M. Z., Multiphase Porous Electrode Theory, [Journal of the Electrochemical Society](https://doi.org/10.1149/2.0171711jes), 2017, 164 (11) E3291-E3310, [arXiv preprint](https://arxiv.org/abs/1702.08432).

## Prerequisites

1.  [Python 3.7](https://www.python.org/) with the following packages installed: `numpy`, `scipy`, `matplotlib`, `pyqt5`, and `h5py`.
2.  [DAE Tools](http://www.daetools.com/) version 1.9.0, which can be [downloaded here](https://sourceforge.net/projects/daetools/files/daetools/1.9.0/).

## Installation

1.  Install the prerequisites above.
2.  Download the [latest release of MPET](https://bitbucket.org/bazantgroup/mpet/downloads/?tab=tags), or clone a copy of this source code repository.
3.  Enter the mpet folder, and use the setup.py script to install the mpet Python package:
    - We recommend using the pip package manager: `pip install .`
    - The legacy approach also works: `python setup.py install`

MPET is also available on [PyPI](https://pypi.org/project/mpet/), the Python Package Index, and can be installed with `pip install mpet`.

## Simulation

1.  Copy the overall system parameters file, `configs/params_system.cfg`, to your working directory.
2.  Copy the material parameter files referred to in the system parameters file (e.g. `configs/params_LFP.cfg` and `configs/params__graphite_1param.cfg`) to the working directory.
3.  Edit `params_system.cfg` to suit the simulation you're trying to run. Be sure to reference a material parameters file for the cathode and optionally one (the same or separate file) for the anode.
4.  Edit the material parameters file(s) serving as the electrode materials.
5.  Run `mpetrun.py`, passing `params_system.cfg` as an argument:
    `mpetrun.py params_system.cfg`

The software will save the simulation output in a time-stamped subdirectory within a directory called `history`. The data contents of the most recent output will also be copied to a directory called `sim_output`. Each output directory should contain:

- the output data (`.mat` file)
- copies of the input parameters files defining the simulation
- a copy of the daetools config parameters (e.g. solver tolerances)
- information about the script used to run the simulation
- information about the simulation (e.g. run time)
- processed, dimensional and nondimensional parameters as
  Python-pickled dictionary objects

## Analysis

1.  Analyze output with `mpetplot.py` (pass output data directory, then plot-type as arguments)
    - e.g., voltage plot: `mpetplot.py sim_output v`
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
2.  Alternatively, convert the output to plain text (csv) format using `mpetplot.py sim_output text` (or replace `sim_output` with any subfolder in the `history` folder). Then analyze using whatever tools you prefer.

If you want to save output to a movie (or figure), add `save` as an extra argument to `mpetplot.py`: `mpetplot.py sim_output cbar save`.

Movie output requires that you have `ffmpeg` or `mencoder` (part of `MPlayer`) installed.

## Troubleshooting

Please use the Issues section of the Bitbucket repository (https://bitbucket.org/bazantgroup/mpet/issues) to file issues and/or bug reports with the software.
