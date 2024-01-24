Parameter space exploration
===========================

The script ``./bin/ensemble.py`` can be used to to automatically generate a set of config files in which one or more parameters can be varied over specified values. The script is executed by running ``./bin/ensemble.py "reference_config_file"``. Within the enseble.py you can indicate which parameter(s) you want to explore over what range. The output will be a number of config files in which all combinations of the given parameter values are used.

Additionally the overarching script ``/bin/mpet_create_runjobs_dashboard.py`` combines the enseble creation, cluster run, and the dashboard. In this script you can idicate the reference config file and the parameter space you want to explore, and the settings of the cluster script you want to use. Next, execute the shell script ``bin/run_mpet_create_run_dashboard.sh``. It will create the config files in the range of inidcated parameters, execute the run cluster script for all created config files, and plot the results of all these models in the dashbaord (which will be started automatically).
