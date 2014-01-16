========
mpet.py
========

1. Install python, numpy, scipy, matplotlib
    - Linux: try using your package manager.
    - Windows: try python(x,y) -- https://code.google.com/p/pythonxy

2. Install DAE Tools
(http://sourceforge.net/projects/daetools/files/1.3.0-beta2/)

3. Copy params_default.cfg to params.cfg

4. Edit params.cfg to suit simulation you're trying to run

5. Run script, passing params.py as an argument:
    - python2 mpet.py params.py

6. Analyze output with plot_data.py (pass output data directory, then
plot-type as arguments).
    - e.g., voltage plot:
        $ python2 plot_data.py sim_output v
    - all options:
      v -- voltage vs filling fraction
      curr -- current vs time
      elytec -- electrolyte concentration (movie)
      elytep -- electrolyte potential (movie)
      csld -- solid concentrations (all, movie)
      phisld -- solid potential (all, used with simSurfCathCond = true, movie)
      cbar -- average solid concentrations (movie)
      cathp -- macroscopic cathode solid phise potential (movie)
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
