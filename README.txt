========
mpet.py
========

1. Install python, numpy, scipy, matplotlib
    - Linux: try using your package manager.
    - Windows: try python(x,y) -- https://code.google.com/p/pythonxy

2. Install DAE Tools
(http://sourceforge.net/projects/daetools/files/1.4.0/)

3. Copy params_default.cfg to params.cfg

4. Edit params.cfg to suit simulation you're trying to run

5. Run script, passing params.py as an argument (try with/without the 2):
    - python[2] mpet.py params.py

The software will save the simulation output data in a folder called
"sim_output" and will also keep a time-stamped copy in a folder called
"history." Each output directory should contain
    - the output data
    - a copy of the input parameters defining the simulation
    - a copy of the daetools config parameters (e.g. solver tolerances)
    - information about the script used to run the simulation
    - information about the simulation (e.g. run time)
    - processed, dimensional and nondimensional parameters as
      Python-pickled dictionary objects

6. Analyze output with plot_data.py (pass output data directory, then
plot-type as arguments).
    - e.g., voltage plot:
        $ python[2] plot_data.py sim_output v
    - all options (full, c, a indicate full cell, cathode, and anode):
      v -- voltage vs filling fraction
      curr -- current vs time
      elytec -- electrolyte concentration (movie)
      elytep -- electrolyte potential (movie)
      csld_{c,a} -- solid concentrations (all, movie, used with solidType_{c,a} not homog)
      phisld_{c,a} -- solid potential (all, movie, used with simSurfCond_{c,a} = true)
      cbar_{full,c,a} -- average solid concentrations (movie)
      bulkp_{c,a} -- macroscopic electrode solid phise potential (movie)
      soc_{c,a} -- electrode state of charge
If you want to save output to a movie (or figure), add "save" (no
quotes) as an extra argument to plot_data.py:
        $ python2 plot_data.py sim_output cbar save
Movie output requires that you have ffmpeg or mencoder (part of
MPlayer) installed.

========
MATLAB Codes (soon to be deprecated)
========
ACR w/ surface wetting and coherency strain
[t,cpcs,ffvec,vvec,disc,part] = acr_mpet_swcs_rev2(dim_crate,part,cpcs0,disc,ffend)

OUTPUTS
t = time
cpcs = results vector
ffvec = filling fraction
vvec = voltage
disc = discretizations
part = particle sizes

INPUTS
dim_crate = dimensional C-rate
part = particle sizes (enter 0 if you dont have one)
cpcs0 = start point (enter 0 if you dont have one, or use last time from previous output)
disc = discretizations (required for use with part vector, enter 0 if not applicable)
ffend = ending filling fraction of discretization (1 for full discharge, 0 for full charge)

Solid solution model with variable particle sizes
[t,cpcs,disc,psd,ffvec,vvec] = pm_2dvarsz_rev1(disc,psd)

OUTPUTS
t = time
cpcs = results vector
disc = discretizations
psd = particle size vector
ffvec = filling fraction
vvec = voltage

INPUTS
disc = discretizations vector (0 if not applicable)
psd = particle sizes (0 if not applicable)

The animate scripts are straightforward.  Just make sure you check the 'fig' and 'output' variables.  Fig denotes the plot variables and output is a boolean (0 or 1) denoting whether or not to record the video.
