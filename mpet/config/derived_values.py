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

        # placeholder for config
        self.config = None

        # store list of available derived values: any method of which
        # the name does not start with _
        # callable is a trick to avoid returning attributes like self.config as well
        self.available_values = [k for k in dir(self) if not k.startswith('_')
                                 and callable(getattr(self, k))]

    def __repr__(self):
        """
        Representation when printing this class:
        print the underlying dict with parameter values
        """
        return dict.__repr__(self.values)

    def __getitem__(self, args):
        """
        Retrieve a derived parameter

        :param tuple args: Tuple with 3 items: Config object, name of parameter
            to retrieve, electrode to retrieve parameter for (None for system values)

        :return: value of derived parameter
        """
        config, item, trode = args

        # set config class-wide for easy access in methods
        self.config = config

        # select general or electrode dict
        if trode is None:
            values = self.values
            func_args = ()
        else:
            values = self.values[trode]
            func_args = (trode, )

        # calculate value if not already stored
        if item not in values:
            try:
                func = getattr(self, item)
            except AttributeError:
                raise UnknownParameterError(f'Unknown parameter: {item}')
            values[item] = func(*func_args)

        return values[item]

    def Damb(self):
        zp = self.config['zp']
        zm = self.config['zm']
        Dp = self.config['Dp']
        Dm = self.config['Dm']
        return ((zp - zm) * Dp * Dm) / (zp * Dp - zm * Dm)

    def tp(self):
        zp = self.config['zp']
        zm = self.config['zm']
        Dp = self.config['Dp']
        Dm = self.config['Dm']
        return zp * Dp / (zp * Dp - zm * Dm)

    def t_ref(self):
        return self.config['L_ref']**2 / self.config['D_ref']

    def curr_ref(self):
        return 3600. / self.config['t_ref']

    def sigma_s_ref(self):
        return self.config['L_ref']**2 * constants.F**2 * constants.c_ref \
            / (self.config['t_ref'] * constants.k * constants.N_A * constants.T_ref)

    def currset(self):
        return self.config['1C_current_density'] * self.config['Crate']  # A/m^2

    def Rser_ref(self):
        return constants.k * constants.T_ref \
            / (constants.e * self.config['curr_ref'] * self.config['1C_current_density'])

    def csmax(self, trode):
        return self.config[trode, 'rho_s'] / constants.N_A

    def cap(self, trode):
        return constants.e * self.config['L'][trode] * (1 - self.config['poros'][trode]) \
            * self.config['P_L'][trode] * self.config[trode, 'rho_s']

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
        if 'a' in self.config['trodes']:
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
        for trode in self.config['trodes']:
            d[trode] = -self.config[trode, 'muR_ref'][0]
        return d
