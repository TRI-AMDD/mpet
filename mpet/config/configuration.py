"""Helper functions for interacting with system parameters.

This module provides functions for various data format exchanges:
 - config files (on disk) <--> dictionaries of parameters (in memory)
 - dictionaries of parameters (in memory) <--> dictionaries of parameters ("pickled" on disk)

It also has various other functions used in processes for things such as generating
distributions from input means and standard deviations.
"""
import configparser
import os
# import pickle

import numpy as np

from mpet.config import schemas
# import mpet.props_am as props_am
from mpet.config.derived_values import DerivedValues
from mpet.exceptions import UnknownParameterError


#: parameter that are define per electrode with a _{electrode} suffix
PARAMS_PER_TRODE = ['Nvol', 'Npart', 'mean', 'stddev', 'cs0', 'simBulkCond', 'sigma_s',
                    'simPartCond', 'G_mean', 'G_stddev', 'L', 'P_L', 'poros', 'BruggExp',
                    'specified_psd']
#: subset of PARAMS_PER_TRODE that is defined for the separator as well
PARAMS_SEPARATOR = ['Nvol', 'L', 'poros', 'BruggExp']
#: parameters that are used with several names. TODO: can we get rid of these?
PARAMS_ALIAS = {'CrateCurr': '1C_current_density', 'n_refTrode': 'n', 'Tabs': 'T',
                'td': 't_ref'}


class Config:
    def __init__(self, paramfile="params.cfg"):
        """
        Hold values from system and electrode configuration files, as well as
        derived values
        """
        # store path to config file
        self.path = os.path.dirname(paramfile)

        # load system parameters file
        self.D_s = ParameterSet(paramfile, 'system', self.path)
        # the anode and separator are optional: only if there are volumes to simulate
        trodes = ['c']
        if self.D_s['Nvol_a'] > 0:
            trodes.append('a')
        have_separator = self.D_s['Nvol_s'] > 0
        self.D_s.params['trodes'] = trodes
        self.D_s.have_separator = have_separator

        # load electrode parameter file(s)
        self.paramfiles = {}

        cathode_paramfile = self.D_s['cathode']
        if not os.path.isabs(cathode_paramfile):
            cathode_paramfile = os.path.join(self.path, cathode_paramfile)
            self.paramfiles['c'] = cathode_paramfile
        self.D_c = ParameterSet(cathode_paramfile, 'electrode', self.path)
        self.D_c.params['trodes'] = trodes
        self.D_c.have_separator = have_separator

        if 'a' in trodes:
            anode_paramfile = self.D_s['anode']
            if not os.path.isabs(anode_paramfile):
                anode_paramfile = os.path.join(self.path, anode_paramfile)
            self.paramfiles['a'] = anode_paramfile
            self.D_a = ParameterSet(anode_paramfile, 'electrode', self.path)
            self.D_a.params['trodes'] = trodes
            self.D_a.have_separator = have_separator
        else:
            self.D_a = None

        # initialize class to calculate and hold derived values
        self.derived_values = DerivedValues()

        # set the random seed
        # TODO: replace by a process_config method that does more processing, e.g.
        # including particle size distribution
        if self.D_s['randomSeed']:
            np.random.seed(self.D_s['seed'])

    def __getitem__(self, items):
        """
        Get the value of a parameter, either a single item to retrieve from the system config,
        or a tuple of (electrode, item), in which case item is read from the config of the given
        electrode (a or c)
        If the value is found in none of the configs, it is assumed to be a derived
        parameter that can be calculated from the config values
        """
        try:
            if isinstance(items, tuple):
                # electrode config
                try:
                    trode, item = items
                except ValueError:
                    raise ValueError(f"Reading from electrode config requires two arguments, but "
                                     f"got {len(items)}")
                # select correct config
                if trode == 'a':
                    assert self.D_a is not None, "Anode parameter requested but " \
                                                 "anode is not simulated"
                    d = self.D_a
                elif trode == 'c':
                    d = self.D_c
                else:
                    raise ValueError(f"Provided electrode must be a or c, got {trode}")
            else:
                # system config
                trode = None
                item = items
                d = self.D_s

            # try to read the parameter from the config
            try:
                value = d[item]
            except UnknownParameterError:
                # not known in config, assume it is a derived value
                # this will raise an UnknownParameterError if still not found
                value = self.derived_values.get(self, item, trode)

        except RecursionError:
            raise Exception(f"Failed to get {items} due to recursion error. "
                            f"Circular parameter dependency?")

        return value

    def __setitem__(self, items, value):
        """
        Set a parameter value, either a single item to store in system config,
        or a tuple of (electrode, item), in which case the item is stored in the config
        of the given electrode (a or c)
        """
        # then process as in __getitem__
        if isinstance(items, tuple):
            try:
                trode, item = items
            except ValueError:
                raise ValueError(f"Setting electrode config value requires two arguments, but "
                                 f"got {len(items)}")
            # select correct config
                if trode == 'a':
                    assert self.D_a is not None, "Anode parameter requested but " \
                                                 "anode is not simulated"
                    d = self.D_a
                elif trode == 'c':
                    d = self.D_c
                else:
                    raise ValueError(f"Provided electrode must be a or c, got {trode}")
        else:
            # system config
            item = items
            d = self.D_s

        d[item] = value


class ParameterSet:
    def __init__(self, paramfile, config_type, path):
        """
        Hold a set of parameters from a single entity (system, electrode)
        """
        assert config_type in ['system', 'electrode'], f"Invalid config type: {config_type}"
        self.path = path
        self.config_type = config_type

        self.have_separator = False
        self.params = {}

        self._load_file(paramfile)

    def _load_file(self, fname):
        """
        Create config from file
        """
        if not os.path.isfile(fname):
            raise Exception(f"Missing config file: {fname}")
        # create config parser
        parser = configparser.ConfigParser()
        parser.optionxform = str
        parser.read(fname)

        # load each section and validate schema
        for section in parser.sections():
            # get the schema
            try:
                config_schema = getattr(schemas, self.config_type)[section]
            except KeyError:
                raise Exception(f"Unknown section '{section}' in {self.fname}")
            # validate
            section_params = config_schema.validate(dict(parser[section].items()))
            # verify there are no duplicate keys
            for key in section_params.keys():
                if key in self.params:
                    raise Exception(f"Duplicate key found: {key}")
            # store the config
            self.params.update(section_params)

    def __repr__(self):
        """
        When printing this class, print the parameters dict
        """
        return dict.__repr__(self.params)

    def __getitem__(self, item):
        """
        Get a parameter
        """
        # if an item is unknown in the dict, try the aliases and per electrode values
        if item not in self.params:
            self.params[item] = self._get_value(item)

        return self.params[item]

    def __setitem__(self, item, value):
        """
        Set a parameter value
        """
        print(f"Setting {item} to {value}")
        # check if the parameter was already set
        # TODO: verify this check is useful, or if overwriting a value is ok
        # note: explicitly check in self.params instead of just self, as otherwise
        # the getter may complain that the parameter does not exist,
        # which is actually what we expect
        if item in self.params:
            raise ValueError(f"Trying to set parameter {item} but it is already defined")

        self.params[item] = value

    def _get_value(self, item):
        """
        Get a value that is defined per electrode/separator, or is used as an alias
        """
        if item in PARAMS_PER_TRODE:
            # create a new dict containg the value per electrode
            d = {}
            # some parameters are also defined for the separator
            trodes = self['trodes'][:]  # make a copy here to avoid adding values to the original
            if item in PARAMS_SEPARATOR and self.have_separator:
                trodes.append('s')
            for trode in trodes:
                # get the value for this electrode/separator and store it
                key = f'{item}_{trode}'
                d[trode] = self[key]
                # delete the original key. TODO: verify this is ok to do
                del self.params[key]
            return d
        elif item in PARAMS_ALIAS:
            return self[PARAMS_ALIAS[item]]
        else:
            raise UnknownParameterError(f"Unknown parameter: {item}")
