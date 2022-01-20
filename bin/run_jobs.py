#!/usr/bin/env python
import argparse
import os
import mpet.main as main

from dask_jobqueue import SLURMCluster, PBSCluster
from dask.distributed import Client, LocalCluster


def create_slurm_cluster(time, nproc, mem, queue, dashboard_port):
    """Create a SLURM cluster for use with dask"""
    cluster = SLURMCluster(cores=nproc,
                           processes=nproc,
                           memory=mem,
                           queue=queue,
                           walltime=time,
                           scheduler_options={'dashboard_address': dashboard_port})
    return cluster


def create_pbs_cluster(time, nproc, mem, queue, dashboard_port):
    """Create a PBS cluster for use with dask"""
    cluster = PBSCluster(cores=nproc,
                         processes=nproc,
                         memory=mem,
                         resource_spec=f'nodes=1:ppn={nproc}',
                         queue=queue,
                         walltime=time,
                         scheduler_options={'dashboard_address': dashboard_port})
    return cluster


#def run_mpet(client, mpet_configs):
#    """Run MPET on each config file present in the mpet_configs folder"""
#    with open(mpet_configs, 'r') as fp:
#        config_files = fp.readlines()

#    folder = os.path.dirname(mpet_configs)

    #commands = [f'. ~/.bashrc && conda activate mpet && cd ~/batteries/output && mpetrun.py '
    #            f'{os.path.join(folder, fname.strip())} -o {fname.strip()} 2>&1 > output.txt'
    #            for fname in config_files]
#    print('Running', commands)
#    futures = client.map(os.system, commands)
    #futures = client.map(main.main('./configs/params_system.cfg'))
#    print('Waiting for MPET to finish')
#    client.gather(futures)
#    print('Done')
def run_mpet(files):
    return main.main(files)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Runs several instances of MPET on a SLURM or'
                                     ' PBS cluster, or local cluster')
    # cluster settings
    parser.add_argument('--time', '-t', required=True,
                        help='Maximum walltime per job (hh:mm:ss format)')
    parser.add_argument('--nproc', '-n', type=int, required=True,
                        help='Number of CPU cores per job')
    parser.add_argument('--mem', '-m', required=True,
                        help='Max memory usage per job')
    parser.add_argument('--queue', '-q', default='default',
                        help='Queue to use (default: %(default)s)')
    parser.add_argument('--dashboard_port', '-d', type=int, default=4096,
                        help='Port for dask dashboard (default: %(default)s)')

    # script settings
    parser.add_argument('--scheduler', '-s', default='slurm', choices=('slurm', 'pbs', 'local'),
                        help='Scheduling system to use (default: %(default)s)')
    parser.add_argument('--min_jobs', type=int, default=1,
                        help='Minimum number of jobs to launch (default: %(default)s)')
    parser.add_argument('--max_jobs', type=int, default=1,
                        help='Maximum number of jobs to launch (default: %(default)s)')
    parser.add_argument('mpet_configs',
                        help='Text file containg the path to each MPET config file to run')

    args = parser.parse_args()
    # split cluster settings from the rest
    main_settings = {}
    for arg in ['scheduler', 'mpet_configs', 'min_jobs', 'max_jobs']:
        main_settings[arg] = getattr(args, arg)
        delattr(args, arg)
    cluster_settings = vars(args)

    # create cluster
    if main_settings['scheduler'] == 'slurm':
        cluster = create_slurm_cluster(**cluster_settings)
    elif main_settings['scheduler'] == 'pbs':
        cluster = create_pbs_cluster(**cluster_settings)
    elif main_settings['scheduler'] == 'local':
        cluster = LocalCluster(n_workers=1, threads_per_worker=1,memory_limit='1GB')

    # Scale Dask cluster automatically based on scheduler activity (only if not local cluster)
    if main_settings['scheduler'] != 'local':
        cluster.adapt(minimum_jobs=main_settings['min_jobs'],
                      maximum_jobs=main_settings['max_jobs'])
    client = Client(cluster)

    mpet_configs = os.path.abspath(main_settings['mpet_configs'])
    with open(mpet_configs, 'r') as fp:
        config_files = fp.readlines()

    folder = os.path.dirname(mpet_configs)
    files = [f'{os.path.join(folder, fname.strip())}' for fname in config_files]

    print('Running', files)
    #main.main(files[0]) # This works
    futures = client.map(run_mpet, files)
    print('Waiting for MPET to finish')
    #input()
    print('going to gather futures')
    client.gather(futures)
    print('Done')

    client.shutdown()
