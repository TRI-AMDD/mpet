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
from mpet.config import constants
from mpet.exceptions import UnknownParameterError


#: parameter that are define per electrode with a _{electrode} suffix
PARAMS_PER_TRODE = ['Nvol', 'Npart', 'mean', 'stddev', 'cs0', 'simBulkCond', 'sigma_s',
                    'simPartCond', 'G_mean', 'G_stddev', 'L', 'P_L', 'poros', 'BruggExp',
                    'specified_psd']
#: subset of PARAMS_PER_TRODE that is defined for the separator as well
PARAMS_SEPARATOR = ['Nvol', 'L', 'poros', 'BruggExp']
#: parameters that are used with several names. TODO: can we get rid of these?
PARAMS_ALIAS = {'CrateCurr': '1C_current_density', 'n_refTrode': 'n', 'Tabs': 'T',
                'td': 't_ref', 'Omga': 'Omega_a', 'Omgb': 'Omega_b', 'Omgc': 'Omega_c'}


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
        self.trodes = ['c']
        if self.D_s['Nvol_a'] > 0:
            self.trodes.append('a')
        self.trodes = self.trodes
        self.have_separator = self.D_s['Nvol_s'] > 0
        self.D_s.params['trodes'] = self.trodes
        self.D_s.have_separator = self.have_separator

        # load electrode parameter file(s)
        self.paramfiles = {}

        cathode_paramfile = self.D_s['cathode']
        if not os.path.isabs(cathode_paramfile):
            cathode_paramfile = os.path.join(self.path, cathode_paramfile)
            self.paramfiles['c'] = cathode_paramfile
        self.D_c = ParameterSet(cathode_paramfile, 'electrode', self.path)
        self.D_c.params['trodes'] = self.trodes
        self.D_c.have_separator = self.have_separator

        if 'a' in self.trodes:
            anode_paramfile = self.D_s['anode']
            if not os.path.isabs(anode_paramfile):
                anode_paramfile = os.path.join(self.path, anode_paramfile)
            self.paramfiles['a'] = anode_paramfile
            self.D_a = ParameterSet(anode_paramfile, 'electrode', self.path)
            self.D_a.params['trodes'] = self.trodes
            self.D_a.have_separator = self.have_separator
        else:
            self.D_a = None

        # initialize class to calculate and hold derived values
        self.derived_values = DerivedValues()

        # set defaults and scale values that should be non-dim
        self.config_processed = False
        self._process_config()

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

            # get alias of item in case it is defined
            try:
                item = PARAMS_ALIAS[item]
            except KeyError:
                # no alias for this item
                pass

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

    def _process_config(self):
        """
        Set default values and process some values
        """
        # TODO: To add in this method: prevDir?
        if self.config_processed:
            raise Exception("The config can be processed only once as values are scaled in-place")
        self.config_processed = True

        # set the random seed (do this before generating any random distributions)
        if self.D_s['randomSeed']:
            np.random.seed(self.D_s['seed'])

        # set default CrateCurr = 1C_current_density
        limtrode = self['limtrode']
        theoretical_1C_current = self[limtrode, 'cap'] / 3600.  # A/m^2
        param = '1C_current_density'
        if self[param] is None:
            # set to theoretical value
            self[param] = theoretical_1C_current

        # non-dimensional scalings
        # TODO: what is an optional parameter needs scaling?
        # first check if parameter is present?
        self['T'] = self['Tabs'] / constants.T_ref
        self['Rser'] = self['Rser'] / self['Rser_ref']
        self['Dp'] = self['Dp'] / self['D_ref']
        self['Dm'] = self['Dm'] / self['D_ref']
        self['c0'] = self['c0'] / constants.c_ref
        self['phi_cathode'] = 0.  # TODO: why is this defined if always 0?
        self['currset'] = self['currset'] / (theoretical_1C_current * self['curr_ref'])
        self['k0_foil'] = self['k0_foil'] / (self['CrateCurr'] * self['curr_ref'])
        self['Rfilm_foil'] = self['Rfilm_foil'] / self['Rser_ref']

        # scalings per electrode
        self['beta'] = {}
        for trode in self.trodes:
            self['L'][trode] = self['L'][trode] / self['L_ref']
            self['beta'][trode] = self[trode, 'csmax'] / constants.c_ref
            self['sigma_s'][trode] = self['sigma_s'][trode] / self['sigma_s_ref']

            kT = constants.k * constants.T_ref
            self[trode, 'lambda'] = self[trode, 'lambda'] / kT
            self[trode, 'B'] = self[trode, 'B'] / (kT * constants.N_A * self[trode, 'cs_ref'])
            for param in ['Omga', 'Omgb', 'Omgc', 'EvdW']:
                value = self[trode, param]
                if value is not None:
                    self[trode, param] = value / kT

        # particle distributions
        # TODO: check for prevDir, load previous values if set
        self._distr_part()

        # Gibss free energy, must be done after distr_part
        self._G()

    def _distr_part(self):
        """
        """
        # intialize dicts in config
        self['psd_num'] = {}
        self['psd_len'] = {}
        self['psd_area'] = {}
        self['psd_vol'] = {}
        self['psr_vol_FracVol'] = {}

        for trode in self.trodes:
            solidType = self[trode, 'type']
            Nvol = self['Nvol'][trode]
            Npart = self['Npart'][trode]

            # check if PSD is specified. If so, it is an ndarray so use np.all
            if not np.all(self['specified_psd'][trode]):
                # If PSD is not specified, make a length-sampled particle size distribution
                # Log-normally distributed
                mean = self['mean'][trode]
                stddev = self['stddev'][trode]
                if np.allclose(stddev, 0., atol=1e-12):
                    raw = mean * np.ones((Nvol, Npart))
                else:
                    var = stddev**2
                    mu = np.log((mean**2) / np.sqrt(var + mean**2))
                    sigma = np.sqrt(np.log(var/(mean**2) + 1))
                    raw = np.random.lognormal(mu, sigma, size=(Nvol, Npart))
            else:
                # use user-defined PSD
                raw = self['specified_psd'][trode]
                if raw.shape != (Nvol, Npart):
                    raise ValueError('Specified particle size distribution discretization '
                                     'of volumes inequal to the one specified in the config file')

            # TODO: need to store raw distribution?

            # For particles with internal profiles, convert psd to
            # integers -- number of steps
            solidDisc = self[trode, 'discretization']
            if solidType in ['ACR']:
                psd_num = np.ceil(raw / solidDisc).astype(int)
                psd_len = solidDisc * psd_num
            elif solidType in ['CHR', 'diffn', 'CHR2', 'diffn2']:
                psd_num = np.ceil(raw / solidDisc).astype(int) + 1
                psd_len = solidDisc * (psd_num - 1)
            # For homogeneous particles (only one 'volume' per particle)
            # ['homog', 'homog_sdn', 'homog2', 'homog2_sdn']
            elif 'homog' in solidType:
                # Each particle is only one volume
                psd_num = np.ones(raw.shape, dtype=np.int)
                # The lengths are given by the original length distr.
                psd_len = raw
            else:
                raise NotImplementedError(f'Unknown solid type: {solidType}')

            # Calculate areas and volumes
            solidShape = self[trode, 'shape']
            if solidShape == 'sphere':
                psd_area = 4 * np.pi * psd_len**2
                psd_vol = (4. / 3) * np.pi * psd_len**3
            elif solidShape == 'C3':
                psd_area = 2 * 1.2263 * psd_len**2
                psd_vol = 1.2263 * psd_len**2 * self[trode, 'thickness']
            elif solidShape == 'cylinder':
                psd_area = 2 * np.pi * psd_len * self[trode, 'thickness']
                psd_vol = np.pi * psd_len**2 * self[trode, 'thickness']
            else:
                raise NotImplementedError(f'Unknown solid shape: {solidShape}')

            # Fraction of individual particle volume compared to total
            # volume of particles _within the simulated electrode
            # volume_
            psd_frac_vol = psd_vol / psd_vol.sum(axis=1, keepdims=True)

            # store values to config
            self['psd_num'][trode] = psd_num
            self['psd_len'][trode] = psd_len
            self['psd_area'][trode] = psd_area
            self['psd_vol'][trode] = psd_vol
            self['psr_vol_FracVol'][trode] = psd_frac_vol

    def _G(self):
        """
        """
        self['G'] = {}
        for trode in self.trodes:
            Nvol = self['Nvol'][trode]
            Npart = self['Npart'][trode]
            mean = self['G_mean'][trode]
            stddev = self['G_stddev'][trode]
            if np.allclose(stddev, 0, atol=1e-12):
                G = mean * np.ones((Nvol, Npart))
            else:
                var = stddev**2
                mu = np.log((mean**2) / np.sqrt(var + mean**2))
                sigma = np.sqrt(np.log(var / (mean**2) + 1))
                G = np.random.lognormal(mu, sigma, size=(Nvol, Npart))

            # scale and store
            self['G'][trode] = G * constants.k * constants.T_ref * self['t_ref'] \
                / (constants.e * constants.F * self[trode, 'csmax'] * self['psd_vol'][trode])


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
        # if an item is unknown in the dict, try the per electrode values
        if item not in self.params:
            self.params[item] = self._get_value(item)

        return self.params[item]

    def __setitem__(self, item, value):
        """
        Set a parameter value
        """
        self.params[item] = value

    def _get_value(self, item):
        """
        Get a value that is defined per electrode/separator
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
        else:
            raise UnknownParameterError(f"Unknown parameter: {item}")