# MPET

1.  Install `python`, `numpy`, `scipy`, `matplotlib`
    - Linux: try using your package manager
    - Windows: try [WinPython](https://winpython.github.io)
2.  Install [DAE Tools](http://sourceforge.net/projects/daetools/files/1.4.0)
3.  Copy the overall system parameters file,
    `configDefaults/params_system.cfg`, to main directory
4.  Copy at least one material parameters file from `configDefaults`
    (e.g. `configDefaults/params_electrode.cfg`) to the main directory
5.  Edit params_system.cfg to suit the simulation you're trying to run. Be
    sure to reference a material parameters file for the cathode and
    optionally one (the same or separate file) for the anode.
6.  Edit the material parameters file(s) serving as the electrode
    materials.
7.  Run script, passing params_system.py as an argument (try with/without the 2):
    - python\[2\] mpet.py params.py

    The software will save the simulation output data in a folder called
    "sim_output" and will also keep a time-stamped copy in a folder called
    "history." Each output directory should contain
    - the output data (.mat file)
    - copies of the input parameters files defining the simulation
    - a copy of the daetools config parameters (e.g. solver tolerances)
    - information about the script used to run the simulation
    - information about the simulation (e.g. run time)
    - processed, dimensional and nondimensional parameters as
      Python-pickled dictionary objects
8.  Analyze output with plot_data.py (pass output data directory, then
    plot-type as arguments).
    - e.g., voltage plot:
        $ python[2] plot_data.py sim_output v
    - other options (full, c, a indicate full cell, cathode, and anode):
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
