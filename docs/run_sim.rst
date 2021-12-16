Run a simulation with MPET
==========================


 * Copy the overall system parameters file, ``configs/params_system.cfg``, to your working directory.
 * Copy the material parameter files referred to in the system parameters file (e.g. ``configs/params_LFP.cfg`` and ``configs/params_graphite_1param.cfg``) to the working directory.
 * Edit ``params_system.cfg`` to suit the simulation you're trying to run. Be sure to reference a material parameters file for the cathode and optionally one (the same or separate file) for the anode.
 * Edit the material parameters file(s) serving as the electrode materials.
 * Run ``mpetrun.py``, passing ``params_system.cfg`` as an argument: ``mpetrun.py params_system.cfg``


The software will save the simulation output in a time-stamped subdirectory within a directory called history. The data contents of the most recent output
will also be copied to a directory called sim_output. Each output directory should contain:

 * the output data (``.mat`` file)
 * an HDF5 file containing the details of the simulation
 * copies of the input parameters files defining the simulation
 * a copy of the daetools config parameters (e.g. solver tolerances)
 * information about the script used to run the simulation
 * information about the simulation (e.g. run time)
 * processed, dimensional and nondimensional parameters as Python-pickled dictionary objects