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
from mpet.config.derived_params import DefinitionsSystem, DefinitionsElectrode


#: parameter that are define per electrode with a _{electrode} suffix
PARAMS_PER_TRODE = ['Nvol', 'Npart', 'mean', 'stddev', 'cs0', 'simBulkCond', 'sigma_s',
                    'simPartCond', 'G_mean', 'G_stddev', 'L', 'P_L', 'poros', 'BruggExp']
#: subset of PARAMS_PER_TRODE that is defined for the separator as well
PARAMS_SEPARATOR = ['Nvol', 'L', 'poros', 'BruggExp']
#: parameters that are used with several names. TODO: can we get rid of these?
PARAMS_ALIAS = {'CrateCurr': '1C_current_density', 'n_refTrode': 'n', 'Tabs': 'T'}


class Config:
    def __init__(self, paramfile="params.cfg"):
        """
        Hold values from system and electrode configuration files
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
        cathode_paramfile = self.D_s['cathode']
        if not os.path.isabs(cathode_paramfile):
            cathode_paramfile = os.path.join(self.path, cathode_paramfile)
        self.D_c = ParameterSet(cathode_paramfile, 'cathode', self.path)
        self.D_c.params['trodes'] = trodes
        self.D_c.have_separator = have_separator

        if 'a' in trodes:
            anode_paramfile = self.D_s['anode']
            if not os.path.isabs(anode_paramfile):
                anode_paramfile = os.path.join(self.path, anode_paramfile)
            self.D_a = ParameterSet(anode_paramfile, 'anode', self.path)
            self.D_a.params['trodes'] = trodes
            self.D_a.have_separator = have_separator
        else:
            self.D_a = None

        # set the random seed
        if self.D_s['randomSeed']:
            np.random.seed(self.D_s['seed'])

    def __getitem__(self, items):
        """
        Get the value of a parameter, either a single item to retrieve from the system config,
        or a tuple of (electrode, item), in which case item is read from the config of the given
        electrode (a or c)
        """
        try:
            if isinstance(items, tuple):
                try:
                    trode, item = items
                except ValueError:
                    raise ValueError(f"Reading from config requires one or two arguments, but "
                                     f"got {len(items)}")
            else:
                trode = None
                item = items

            if trode is None:
                # get parameter from system config
                return self.D_s[item]
            else:
                # get parameter from electrode config
                if trode == 'a':
                    assert self.D_a is not None, "Anode parameter requested but " \
                                                 "anode is not simulated"
                    return self.D_a[item]
                elif trode == 'c':
                    return self.D_c[item]
                else:
                    raise ValueError(f"Provided electrode must be a or c, got {trode}")
        except RecursionError:
            raise Exception(f"Failed to get {items} due to recursion error. "
                            f"Circular parameter dependency?")


class ParameterSet:
    def __init__(self, paramfile, entity, path):
        """
        Hold a set of parameters from a single entity (system, cathode, anode)
        """
        assert entity in ['system', 'cathode', 'anode'], f"Invalid entity: {entity}"
        self.path = path

        if entity == 'system':
            self.definitions = DefinitionsSystem()
            self.config_type = 'system'
        else:
            # cathode or anode
            self.definitions = DefinitionsElectrode(entity)
            # config type is the same for cathode/anode
            self.config_type = 'electrode'

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
        Get a parameter, either from known parameters or a derived value
        """
        # if an item is unknown in the dict, assume it needs to be calculated
        if item not in self.params:
            value = self._get_value(item)
            if value is None:
                raise Exception(f"Unknown parameter: {item}")
            # store the calculated value
            self.params[item] = value

        return self.params[item]

    def _get_value(self, item):
        """
        Calculate a value that is not part of config file
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
            # assume this is a derived parameter
            # return None when not found, which is further
            # handled in __getitem__
            # TODO: Some values may actually be None, perhaps handle through exceptions instead
            return self.definitions.get(item, self)
        # TODO: properly handle prevDir with the new schema format
        # elif item == 'prevDir':
        #     value = self.parser.get(item)
        #     if value != 'false' and not os.path.isabs(value):
        #         return os.path.normpath(os.path.join(self.path, value))
