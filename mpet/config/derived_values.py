from string import Formatter as StringFormatter

import numpy as np

from mpet import props_elyte, props_am
from mpet.exceptions import UnknownParameterError
from mpet.config import constants


class DerivedValues:
    def __init__(self):
        """
        DerivedValues holds the functions and equations required
        to calculate a parameter that is derived from other parameters.
        Results are cached, so each parameter is calculated only once
        """
        # to keep track of values that were already calculated
        # initialize with empty dicts for electrodes
        self.values = {'c': {}, 'a': {}}

        self.config = None

        # Equations defining derived parameters
        # any parameter must be enclosed in curly brackets
        # electrode-specific parameters must have a _tr suffix
        # electrode selection from system config must have [tr] suffix
        # constants can be accessed with constants.<constants name>, without curly brackets
        self.equations = {'Damb': '(({zp} - {zm}) * {Dp} * {Dm}) / ({zp} * {Dp} - {zm} * {Dm})',
                          'tp': '{zp} * {Dp} / ({zp} * {Dp} - {zm} * {Dm})',
                          't_ref': '{L_ref}**2 / {D_ref}',
                          'curr_ref': '3600. / {t_ref}',
                          'sigma_s_ref': '{L_ref}**2 * constants.F**2 * constants.c_ref '
                                         '/ ({t_ref} * constants.k * constants.N_A '
                                         '* constants.T_ref)',
                          'currset': '{CrateCurr} * {Crate}',  # A/m^2
                          'Rser_ref': 'constants.k * constants.T_ref / (constants.e * {curr_ref} '
                                      '* {CrateCurr})',
                          'csmax': '{rho_s_tr} / constants.N_A',
                          'cap': 'constants.e * {L[tr]} * (1 - {poros[tr]}) '
                                 ' * {P_L[tr]} * {rho_s_tr}'}

    def __repr__(self):
        """
        Representation when printing this class:
        print the underlying dict with parameter values
        """
        return dict.__repr__(self.values)

    def get(self, config, item, trode=None):
        """
        Retrieve a derived parameter

        :param Config config: Global configuration object
        :param str item: Name of parameter
        :param str trode: Electrode to retrieve parameter for (None for system values)
        """
        # set config class-wide for easy access in methods
        self.config = config

        # select general or electrode dict
        if trode is None:
            values = self.values
            args = ()
        else:
            values = self.values[trode]
            args = (trode, )

        # calculate value if not already stored
        if item not in values:
            if item in self.equations:
                values[item] = self._process_equation(self.equations[item], trode)
            else:
                # get the method to calculate the value
                try:
                    func = getattr(self, item)
                except AttributeError:
                    raise UnknownParameterError(f'Unknown parameter: {item}')
                values[item] = func(*args)

        return values[item]

    def _process_equation(self, equation, trode=None):
        """
        Calculate a parameter that is defined in an equation defined as string

        :param str equation: Equation to evaluate
        :param str trode: electrode, None for system values
        """
        # get parameters that are specified in equation, i.e. everything in curly brackets
        params = [item[1] for item in StringFormatter().parse(equation) if item[1] is not None]
        formatter = {}
        # add tr key to formatter in case param[tr] is used in the string formatting
        if trode is not None:
            formatter['tr'] = trode
        have_electrode_select = False
        for param in params:
            # suffixes for electrode parameters
            electrode_suffix = '_tr'
            electrode_select = '[tr]'
            if param.endswith(electrode_suffix):
                # electrode parameter
                assert trode is not None, 'Equation requires electrode specification'
                # remove trailing _tr from param. Don't use .split in case _tr occurs
                # somewhere in the middle as well
                param_config = param[:-len(electrode_suffix)]
                formatter[param] = self.config[trode, param_config]
            elif param.endswith(electrode_select):
                # electrode selected from system parameter
                assert trode is not None, 'Equation requires electrode specification'
                have_electrode_select = True
                # remove trailing [tr]
                param_config = param[:-len(electrode_select)]
                # put the full parameter dict in the formatter, as the []
                # operator still works in string formatting
                formatter[param_config] = self.config[param_config]
            else:
                # system parameter
                formatter[param] = self.config[param]
        # if there is electrode selection, replace the tr placeholder
        # by the electrode to use
        if have_electrode_select:
            equation = equation.replace('[tr]', f'[{trode}]')
        return eval(equation.format(**formatter))

    def numsegments(self):
        """
        Number of segments
        """
        return len(self.config['segments'])

    def L_ref(self):
        """
        reference L
        """
        return self.config['L']['c']

    def D_ref(self):
        """
        reference D
        """
        if self.config['elyteModelType'] == 'dilute':
            return self.config['Damb']
        else:
            return getattr(props_elyte, self.config['SMset'])()[-1]

    def z(self):
        """
        z
        """
        if 'a' in self.config.trodes:
            return self.config['c', 'cap'] / self.config['a', 'cap']
        else:
            # flat plate anode with assumed infinite supply of metal
            return 0.

    def limtrode(self):
        """
        limtrode
        """
        if self.config['z'] < 1:
            return 'c'
        else:
            return 'a'

    def cs_ref(self, trode):
        """
        reference cs
        """
        if self.config[trode, 'type'] in constants.one_var_types:
            prefac = 1
        elif self.config[trode, 'type'] in constants.two_var_types:
            prefac = .5
        return prefac * self.config[trode, 'csmax']

    def muR_ref(self, trode):
        """
        reference muR
        """
        muRfunc = props_am.muRfuncs(self.config, trode).muRfunc
        cs0bar = self.config['cs0'][trode]
        cs0 = np.array([cs0bar])

        solidType = self.config[trode, 'type']
        if solidType in constants.two_var_types:
            muR_ref = -muRfunc((cs0, cs0), (cs0bar, cs0bar), 0.)[0][0]
        elif solidType in constants.one_var_types:
            muR_ref = -muRfunc(cs0, cs0bar, 0.)[0]
        else:
            raise ValueError(f'Unknown solid type: {solidType}')
        return muR_ref

    def phiRef(self):
        """
        reference phi
        """
        d = {}
        for trode in self.config.trodes:
            d[trode] = -self.config[trode, 'muR_ref'][0]
        return d
