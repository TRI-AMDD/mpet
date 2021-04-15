"""Helper functions for interacting with system parameters.

This module provides functions for various data format exchanges:
 - config files (on disk) <--> dictionaries of parameters (in memory)
 - dictionaries of parameters (in memory) <--> dictionaries of parameters ('pickled' on disk)

It also has various other functions used in processes for things such as generating
distributions from input means and standard deviations.
"""
import configparser
import os
import pickle

import numpy as np

from mpet.config import schemas
from mpet.config.constants import PARAMS_ALIAS, PARAMS_PER_TRODE, PARAMS_SEPARATOR

from mpet.config.derived_values import DerivedValues
from mpet.config import constants
from mpet.exceptions import UnknownParameterError
# temporary:
from mpet.io_utils import size2regsln


class Config:
    def __init__(self, paramfile='params.cfg', from_dicts=False):
        """
        Hold values from system and electrode configuration files, as well as
        derived values
        """
        # initialize class to calculate and hold derived values
        self.derived_values = DerivedValues()
        self.params_per_particle = []
        if from_dicts:
            # read existing dictionaries instead of parameter file
            # paramfile is now folder with input dicts
            self.path = os.path.normpath(paramfile)
            # set which electrodes there are based on which dict files exist
            self.trodes = ['c']
            if os.path.isfile(os.path.join(self.path, 'input_dict_anode.p')):
                self.trodes.append('a')
            # create empty parametersets and derived values
            self.D_s = ParameterSet(None, 'system', self.path)
            self.D_c = ParameterSet(None, 'electrode', self.path)
            if 'a' in self.trodes:
                self.D_a = ParameterSet(None, 'electrode', self.path)
            else:
                self.D_a = None
            # now populate the dicts
            self.read(self.path, full=True)
            if 's' in self.D_s['Nvol']:
                self.have_separator = True
            else:
                self.have_separator = False
            self.D_s.have_separator = self.have_separator
            self.D_c.have_separator = self.have_separator
            self.D_s.trodes = self.trodes
            self.D_c.trodes = self.trodes
            if self.D_a is not None:
                self.D_a.have_separator = self.have_separator
                self.D_a.trodes = self.trodes
            self.params_per_particle = list(constants.PARAMS_PARTICLE.keys())
        else:
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

            # set defaults and scale values that should be non-dim
            self.config_processed = False
            # either process the config, or read already processed config from disk
            prevDir = self['prevDir']
            if prevDir and prevDir != 'false':
                if not os.path.isabs(prevDir):
                    # assume it is relative to the input parameter files
                    self['prevDir'] = os.path.normpath(os.path.join(self.path, prevDir))
                self._process_config(self['prevDir'])
            else:
                self._process_config()

    @classmethod
    def from_dicts(cls, path):
        """
        Create a config instance from a set of dictionaries on disk, instead
        of from config files
        """
        return cls(path, from_dicts=True)

    def _select_config(self, items):
        """
        Select system or electrode config based on reqested parameter(s)
        """
        if isinstance(items, tuple):
            # electrode config
            try:
                trode, item = items
            except ValueError:
                raise ValueError(f'Reading from electrode config requires two arguments, but '
                                 f'got {len(items)}')
            # select correct config
            if trode == 'a':
                assert self.D_a is not None, 'Anode parameter requested but ' \
                                             'anode is not simulated'
                d = self.D_a
            elif trode == 'c':
                d = self.D_c
            else:
                raise ValueError(f'Provided electrode must be a or c, got {trode}')
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

        # select individual particle settings if requested
        if item in self.params_per_particle:
            assert trode is not None, 'Particle-specific parameter requested but ' \
                                      'electrode not specified'
            d = self[trode, 'indvPart']
        return d, item, trode

    def write(self, folder=None):
        """
        Write config to disk
        """
        filenamebase = 'input_dict'
        if folder:
            filenamebase = os.path.join(folder, filenamebase)

        # system, derived values, cathode, optionally anode
        dicts = {'system': self.D_s.params, 'derived_values': self.derived_values.values,
                 'cathode': self.D_c.params}
        if 'a' in self.trodes:
            dicts['anode'] = self.D_a.params

        for section, d in dicts.items():
            with open(f'{filenamebase}_{section}.p', 'wb') as f:
                pickle.dump(d, f)

    def read(self, folder=None, full=False):
        """
        Read previously processed config from disk
        """
        filenamebase = 'input_dict'
        if folder:
            filenamebase = os.path.join(folder, filenamebase)

        # system, derived values, cathode, optionally anode
        sections = ['system', 'derived_values', 'cathode']
        if 'a' in self.trodes:
            sections.append('anode')

        for section in sections:
            with open(f'{filenamebase}_{section}.p', 'rb') as f:
                try:
                    d = pickle.load(f)
                except UnicodeDecodeError:
                    d = pickle.load(f, encoding='latin1')

            if full:
                # update all config
                if section == 'system':
                    self.D_s.params = d
                elif section == 'derived_values':
                    self.derived_values.values = d
                elif section == 'cathode':
                    self.D_c.params = d
                elif section == 'anode':
                    self.D_a.params = d
                else:
                    raise Exception(f'Unknown section: {section}')
            else:
                # only update generated distributions
                if section == 'system':
                    for key in ['psd_num', 'psd_len', 'psd_area', 'psd_vol',
                                'psd_vol_FracVol', 'G']:
                        self[key] = d[key]
                elif section in ['anode', 'cathode']:
                    trode = section[0]
                    self[trode, 'indvPart'] = d['indvPart']

        # make sure to set the numpy random seed, as is usually done in process_config
        if self['randomSeed']:
            np.random.seed(self.D_s['seed'])

        self.params_per_particle = list(constants.PARAMS_PARTICLE.keys())
        self.config_processed = True

    def __getitem__(self, items):
        """
        Get the value of a parameter, either a single item to retrieve from the system config,
        or a tuple of (electrode, item), in which case item is read from the config of the given
        electrode (a or c)
        If the value is found in none of the configs, it is assumed to be a derived
        parameter that can be calculated from the config values
        """
        # select config
        d, item, trode = self._select_config(items)

        try:
            value = d[item]
        except UnknownParameterError:
            # not known in config, assume it is a derived value
            # this will raise an UnknownParameterError if still not found
            value = self.derived_values.get(self, item, trode)

        return value

    def __setitem__(self, items, value):
        """
        Set a parameter value, either a single item to store in system config,
        or a tuple of (electrode, item), in which case the item is stored in the config
        of the given electrode (a or c)
        """
        # select config, ignore returned trode value as we don't need it
        d, item, _ = self._select_config(items)
        d[item] = value

    def __delitem__(self, items):
        """
        Delete a parameter value, either a single item to delete from system config,
        or a tuple of (electrode, item), in which case the item is deleted from the config
        of the given electrode (a or c)
        """
        # select config, ignore returned trode value as we don't need it
        d, item, _ = self._select_config(items)

        del d[item]

    def _process_config(self, prevDir=None):
        """
        Set default values and process some values
        """
        if self.config_processed:
            raise Exception('The config can be processed only once as values are scaled in-place')
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
        # first check if parameter is present?
        kT = constants.k * constants.T_ref
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

            self[trode, 'lambda'] = self[trode, 'lambda'] / kT
            self[trode, 'B'] = self[trode, 'B'] / (kT * constants.N_A * self[trode, 'cs_ref'])
            for param in ['Omga', 'Omgb', 'Omgc', 'EvdW']:
                value = self[trode, param]
                if value is not None:
                    self[trode, param] = value / kT
        # scalings on separator
        if self.have_separator:
            self['L']['s'] /= self['L_ref']

        # scaling/addition of macroscopic input information
        Vref = self['phiRef']['c']
        if 'a' in self.trodes:
            Vref -= self['phiRef']['a']
        factor = constants.e / kT
        self['Vset'] = -(factor * self['Vset'] + Vref)
        self['phimin'] = -(factor * self['Vmax'] + Vref)
        self['phimax'] = -(factor * self['Vmin'] + Vref)

        # Scaling of current and voltage segments
        segments = []
        if self['profileType'] == 'CCsegments':
            for i in range(len(self['segments'])):
                segments.append((self["segments"][i][0]*self["CrateCurr"]
                                / theoretical_1C_current/self['curr_ref'],
                                self["segments"][i][1]*60/self['t_ref']))
        elif self['profileType'] == 'CVsegments':
            for i in range(len(self['segments'])):
                segments.append((-((constants.e/kT)*self['segments'][i][0]+Vref),
                                self['segments'][i][1]*60/self['t_ref']))

        # Current or voltage segments profiles
        segments_tvec = np.zeros(2 * self['numsegments'] + 1)
        segments_setvec = np.zeros(2 * self['numsegments'] + 1)
        if self['profileType'] == 'CCsegments':
            segments_setvec[0] = 0.
        elif self['profileType'] == 'CVsegments':
            segments_setvec[0] = -(kT / constants.e) * Vref
        tPrev = 0.
        for segIndx in range(self['numsegments']):
            tNext = tPrev + self['tramp']
            segments_tvec[2*segIndx+1] = tNext
            tPrev = tNext
            # Factor of 60 here to convert to s
            tNext = tPrev + (self['segments'][segIndx][1] * 60 - self["tramp"])
            segments_tvec[2*segIndx+2] = tNext
            tPrev = tNext
            setNext = self['segments'][segIndx][0]
            segments_setvec[2*segIndx+1] = setNext
            segments_setvec[2*segIndx+2] = setNext
        segments_tvec /= self['t_ref']
        if self['profileType'] == 'CCsegments':
            segments_setvec /= self['curr_ref']
        elif self['profileType'] == 'CVsegments':
            segments_setvec = -((constants.e/kT)*segments_setvec + Vref)
        if 'segments' in self['profileType']:
            self['tend'] = segments_tvec[-1]
            # Pad the last segment so no extrapolation occurs
            # TODO: uncomment when fixed in dev branch (updated value is not used, breaking test)
            # segments_tvec[-1] *= 1.01
        else:
            self['tend'] /= self['t_ref']
        self['segments'] = segments
        self['segments_tvec'] = segments_tvec
        self['segments_setvec'] = segments_setvec
        if self['profileType'] == 'CC' and not np.allclose(self['currset'], 0., atol=1e-12):
            self['tend'] = np.abs(self['capFrac'] / self['currset'])

        if prevDir:
            # load particle distrubtions etc. from previous run
            self.read(prevDir, full=False)
            # Set params per particle as would be done in _invdPart when generating distributions
            self.params_per_particle = list(constants.PARAMS_PARTICLE.keys())
        else:
            # particle distributions
            self._distr_part()
            # Gibss free energy, must be done after distr_part
            self._G()
            # Electrode parameters that depend on invidividual particle
            self._indvPart()

    def _distr_part(self):
        """
        """
        # intialize dicts in config
        self['psd_num'] = {}
        self['psd_len'] = {}
        self['psd_area'] = {}
        self['psd_vol'] = {}
        self['psd_vol_FracVol'] = {}

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
            elif solidType in ['homog', 'homog_sdn', 'homog2', 'homog2_sdn']:
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
            self['psd_vol_FracVol'][trode] = psd_frac_vol

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

    def _indvPart(self):
        """
        """
        for trode in self.trodes:
            Nvol = self['Nvol'][trode]
            Npart = self['Npart'][trode]
            self[trode, 'indvPart'] = {}

            # intialize parameters
            # TODO: verify dtypes

            for param, dtype in constants.PARAMS_PARTICLE.items():
                self[trode, 'indvPart'][param] = np.empty((Nvol, Npart), dtype=dtype)

            # reference scales per trode
            cs_ref_part = constants.N_A * self[trode, 'cs_ref']  # part/m^3

            # calculate the values for each particle in each volume
            for i in range(Nvol):
                for j in range(Npart):
                    # This specific particle dimensions
                    self[trode, 'indvPart']['N'][i, j] = self['psd_num'][trode][i,j]
                    plen = self['psd_len'][trode][i,j]
                    parea = self['psd_area'][trode][i,j]
                    pvol = self['psd_vol'][trode][i,j]
                    # Define a few reference scales
                    F_s_ref = plen * cs_ref_part / self['t_ref']  # part/(m^2 s)
                    i_s_ref = constants.e * F_s_ref  # A/m^2
                    kappa_ref = constants.k * constants.T_ref * cs_ref_part * plen**2  # J/m
                    gamma_S_ref = kappa_ref / plen  # J/m^2
                    # non-dimensional quantities
                    kappa = self[trode, 'indvPart']['kappa'][i, j] \
                        = self[trode, 'kappa'] / kappa_ref
                    nd_dgammadc = self[trode, 'dgammadc'] * cs_ref_part / gamma_S_ref
                    self[trode, 'indvPart']['beta_s'][i, j] = nd_dgammadc / kappa
                    self[trode, 'indvPart']['D'][i, j] = self[trode, 'D'] * self['t_ref'] / plen**2
                    self[trode, 'indvPart']['E_D'][i, j] = self[trode, 'E_D'] \
                        / (constants.k * constants.N_A * constants.T_ref)
                    self[trode, 'indvPart']['k0'][i, j] = self[trode, 'k0'] \
                        / (constants.e * F_s_ref)
                    self[trode, 'indvPart']['E_A'][i, j] = self[trode, 'E_A'] \
                        / (constants.k * constants.N_A * constants.T_ref)
                    self[trode, 'indvPart']['Rfilm'][i, j] = self[trode, 'Rfilm'] \
                        / (constants.k * constants.T_ref / (constants.e * i_s_ref))
                    self[trode, 'indvPart']['delta_L'][i, j] = (parea * plen) / pvol
                    # If we're using the model that varies Omg_a with particle size,
                    # overwrite its value for each particle
                    if self[trode, 'type'] in ['homog_sdn', 'homog2_sdn']:
                        self[trode, 'indvPart']['Omega_a'][i, j] = size2regsln(plen)
                    else:
                        # just use global value
                        self[trode, 'indvPart']['Omega_a'][i, j] = self[trode, 'Omega_a']

        # store which items are defined per particle, so in the future they are retrieved
        # per particle instead of from the values per electrode
        self.params_per_particle = list(constants.PARAMS_PARTICLE.keys())


class ParameterSet:
    def __init__(self, paramfile, config_type, path):
        """
        Hold a set of parameters from a single entity (system, electrode)
        """
        assert config_type in ['system', 'electrode'], f'Invalid config type: {config_type}'
        self.path = path
        self.config_type = config_type

        self.have_separator = False
        self.params = {}

        if paramfile is not None:
            self._load_file(paramfile)

    def _load_file(self, fname):
        """
        Create config from file
        """
        if not os.path.isfile(fname):
            raise Exception(f'Missing config file: {fname}')
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
                raise Exception(f'Unknown section "{section}" in {self.fname}')
            # validate
            section_params = config_schema.validate(dict(parser[section].items()))
            # verify there are no duplicate keys
            for key in section_params.keys():
                if key in self.params:
                    raise Exception(f'Duplicate key found: {key}')
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
        if item in self.params:
            return self.params[item]
        elif item in PARAMS_PER_TRODE:
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
                del self.params[key]
            self.params[item] = d
            return d
        else:
            raise UnknownParameterError(f'Unknown parameter: {item}')

    def __setitem__(self, item, value):
        """
        Set a parameter value
        """
        self.params[item] = value

    def __delitem__(self, item):
        """
        Remove a parameter value
        """
        del self.params[item]
