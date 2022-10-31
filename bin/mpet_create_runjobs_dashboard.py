from create_ensemble import create_ensemble
# import bin.run_jobs
# import mpet_plot_app
import configparser
import sys


if __name__ == '__main__':
    # Read in file
    if len(sys.argv) < 2:
        print("need the config file [python create_enamble.py <baseconfig>]")
        exit(1)
    cfg = configparser.ConfigParser()
    cfg.optionxform = str
    cfg.read(sys.argv[1])
    create_ensemble(cfg)
