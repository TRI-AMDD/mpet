from create_ensemble import create_ensemble
from run_jobs import create_slurm_cluster, create_pbs_cluster, create_local_cluster, run_mpet
# import mpet_plot_app
import os
import subprocess
import shutil
from dask.distributed import Client


# --------------------------------------------------------------
# Fill all these out
# --------------------------------------------------------------
''' Ensemble settings/ '''
# The default config file that you want to adjust in the ensemble
# (path relative to main folder)
cff = 'configs/params_system.cfg'
# Values that need ensemble
ensemble = [
    [("Sim Params","Nvol_c"), ["4", "5"]]
]


class ClusterSettings():
    '''Cluster settings.'''
    # time = 00:00:00  # Max walltime per job (hh:mm:ss format). Argument is not used with a local.
    nproc = 1  # type=int, Number of CPU cores per job. Argument is not used with a local cluster.
    mem = 1  # Max memory usage per job. For alocal cluster it sets the memory limit per worker.
    queue = 1  # Queue to use. Argument is not used with a local cluster.
    dashboard_port = 4096  # Port for dask dashboard


class MainSettings():
    '''Cluster script settings.'''
    scheduler = 'local'  # choices=('slurm', 'pbs', 'local'); Scheduling system to use
    min_jobs = 1  # int; Minimum number of jobs to launch. Argument is not used with local cluster.
    max_jobs = 1  # int; Maximum number of jobs to launch. Argument is not used with local cluster.
    # Text file containg the path to each MPET config file to run;
    # 'parallel_configs.txt' is what create_ensemble saves to automatically (so don't change it)
    mpet_configs = 'ensemble_parallel_configs.txt'


# --------------------------------------------------------------
'''
End of parameters to fill out.
To run: execute run_mpet_create_run_dashboard.sh from the terminal;
. ./run_mpet_create_run_dashboard.sh
To see the resulting dashboard open the http link presented in the terminal:
" Dash is running on http://127.0.0.1:8050/ "
'''
# --------------------------------------------------------------


# helpers for ensemble, do not change
keys = [vals[0] for vals in ensemble]
val = [vals[1] for vals in ensemble]


def call_run_cluster(output_folder):
    # split cluster settings from the rest
    args = ClusterSettings()
    main_settings = MainSettings()
    '''
    for arg in ["scheduler", "mpet_configs", "min_jobs", "max_jobs"]:
        main_settings[arg] = getattr(args, arg)
        delattr(args.__class__, arg)'''
    cluster_settings = vars(args)

    # create cluster
    if main_settings.scheduler == 'slurm':
        cluster = create_slurm_cluster(**cluster_settings)
    elif main_settings.scheduler == 'pbs':
        cluster = create_pbs_cluster(**cluster_settings)
    elif main_settings.scheduler == 'local':
        cluster = create_local_cluster(args.mem, args.dashboard_port)

    # Scale Dask cluster automatically based on scheduler activity (only if not local cluster)
    if main_settings.scheduler != 'local':
        cluster.adapt(minimum_jobs=main_settings['min_jobs'],
                      maximum_jobs=main_settings['max_jobs'])
    client = Client(cluster)

    run_mpet(client, output_folder, os.path.abspath(main_settings.mpet_configs))
    return


if __name__ == '__main__':
    # Read in config file
    create_ensemble(cff, keys, val)

    # Define output folder
    # Store output in folder this script was called from
    output_folder = './runjobs_dashboard'
    # remove sim output if it already exists to only keep newest output
    if os.path.exists(os.path.join(output_folder, 'sim_output')):
        shutil.rmtree(os.path.join(output_folder, 'sim_output'))
    # create output folder if it does not exist yet
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
    call_run_cluster(output_folder)
    subprocess.call(["python", "./bin/mpet_plot_app.py", "-d",
                     str(os.path.join(output_folder, 'sim_output'))])
