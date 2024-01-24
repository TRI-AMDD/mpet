Running multiple simulations on a cluster
=========================================

If you have many simulations you want to run, you can use ``bin/run_jobs.py`` to run them efficiently on a cluster using `Dask <https://dask.org>`_, either locally or on a slurm or PBS cluster. Using the parallel running option requires the following package to be installed: ``dask-jobqueue``.

1. Follow steps 1-4 from :doc:`run_sim` for each of the simulations you want to run. Then create a text file in your working directory containing the system parameter files for your simulations. This text file should contain the file names of each of the system parameter configuration files for which you want to run a simulation. For example, if you have all your parameter files saved in the ``configs`` directory, create: ``configs/parallel_configs.txt``, which contains the lines::

    params_system.cfg
    params_system_XX.cfg
    params_system_YY.cfg
    etc.

2. Run multiple simulations on a cluster using ``run_jobs.py``. The simplest way to run it, is to run the script on the login node. Pass the text file containing the system parameter files (e.g. ``configs/parallel_configs.txt``) and the cluster arguments:

    - ``-s``: scheduler type. Options: ``slurm``, ``pbs``, and ``local``. Default is ``slurm``.
    - ``-t``: Maximum walltime per job (hh:mm:ss format). Argument is not used with a local cluster.
    - ``-n``: Number of CPU cores and instances of MPET per job. Argument is not used with a local cluster.
    - ``-m``: Max memory usage per job (e.g. 2GB). When using a local cluster it sets the memory limit per worker process.
    - ``-q``: Queue to use. Argument is not used with a local cluster.
    - ``-d``: Port for Dask dashboard (default 4096).
    - ``--min_jobs``: Minimum number of jobs to launch. Default = 1. Argument is not used with a local cluster.
    - ``--max_jobs``: Maximum number of jobs to launch. Default = 1. Argument is not used with a local cluster.
3. The simulation output is the same as described above. For each simulation a separate output folder is created in the ``history`` folder.
