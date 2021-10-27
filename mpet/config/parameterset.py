import os
import configparser

from mpet.config import schemas
from mpet.config.constants import PARAMS_PER_TRODE, PARAMS_SEPARATOR
from mpet.exceptions import UnknownParameterError


class ParameterSet:
    def __init__(self, paramfile, config_type, path):
        """
        Hold a set of parameters for a single entity (system, one electrode).

        :param str paramfile: Full path to .cfg file on disk
        :param str config_type: "system" or "electrode"
        :param str path: Folder containing the .cfg file

        .. make sure to document methods related to [] operator
        .. automethod:: __getitem__
        .. automethod:: __setitem__
        .. automethod:: __delitem__
        """
        assert config_type in ['system', 'electrode'], f'Invalid config type: {config_type}'
        self.path = path
        self.config_type = config_type

        self.params = {}

        if paramfile is not None:
            self._load_file(paramfile)

    def _load_file(self, fname):
        """
        Create config from file.

        :param str fname: full path to .cfg file
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
                raise Exception(f'Unknown section "{section}" in {fname}')
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
        When printing this class, print the parameters dict.
        """
        return dict.__repr__(self.params)

    def __getitem__(self, item):
        """
        Get a parameter using the ``[]`` operator.
        See ``Config.__getitem__`` for example usage of the ``[]`` operator.

        :param str item: Name of the parameter

        :return: value of the parameter
        """
        if item in self.params:
            return self.params[item]
        elif item in PARAMS_PER_TRODE:
            # this is only valid for electrode config
            assert self.config_type == 'system', 'Requested electrode parameter from system config'
            # create a new dict containg the value per electrode
            d = {}
            # some parameters are also defined for the separator
            trodes = self['trodes'][:]  # make a copy here to avoid adding values to the original
            if item in PARAMS_SEPARATOR and self['have_separator']:
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
        Set a parameter value.

        :param str item: Name of the parameter
        :param value: Value of the parameter
        """
        self.params[item] = value

    def __delitem__(self, item):
        """
        Remove a parameter value.

        :param str item: Name of the parameter
        """
        del self.params[item]
