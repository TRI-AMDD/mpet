import numpy as np

from mpet import props_elyte, props_am
from mpet.exceptions import UnknownParameterError
from mpet.config import constants


class DerivedValues:
    def __init__(self):
        """
        DerivedValues holds the functions and equations required
        to calculate a parameter that is derived from other parameters.
        Results are cached, so each parameter is calculated only once.
        This class is meant to be used from within :class:`mpet.config.configuration.Config`.

        Like :class:`mpet.config.configuration.Config`, the ``[]`` operator
        can be used to access values.

        Each method of this class that does not start with an underscore
        is a derived value that can be calculated. If it has a ``trode``
        argument, the parameter is electrode-specific and the name
        of the electrode should be specified (either 'a' or 'c').

        A list of the available parameters is stored in ``self.available_values``.

        Where applicable, all parameters are scaled to their non-dimensional values.

        .. make sure to document methods related to [] operator
        .. automethod:: __getitem__
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
        Retrieve a derived parameter using the ``[]`` operator.

        :param tuple args: Tuple with 2 or 3 items: Config object, name of parameter
            to retrieve, optionally electrode to retrieve parameter for (system
            config is selected if no electrode is provided or value is set to None)

        :return: value of derived parameter

        Example usage:

        Initialize config and derived values objects:
        Note that ``DerivedValues`` is meant to be used from within the
        :class:`mpet.config.configuration.Config` class. Passing the entire
        ``Config`` object is then achieved by passing ``self``, instead of ``config``
        to ``DerivedValues``.

        >>> config = Config('/path/to/params_system.cfg')
        >>> dv = DerivedValues()

        Parameter that does not need electrode, ``t_ref``:

        >>> dv[config, 't_ref', None]
        7.792736808950305

        Parameter that does need electrode:

        >>> dv[config, 'csmax', 'c']
        22904.35071404849

        If a parameter that needs and electrode is accessed without electrode,
        or vice-versa, an error is raised:

        >>> dv[config, 'csmax', None]
        TypeError: csmax() missing 1 required positional argument: 'trode'

        >>> dv[config, 't_ref', 'c']
        TypeError: t_ref() takes 1 positional argument but 2 were given
        """

        if len(args) == 2:
            config, item = args
            trode = None
        elif len(args) == 3:
            config, item, trode = args
        else:
            raise ValueError(f"There should be 2 or 3 arguments, got {len(args)}")

        # set config class-wide for easy access in methods
        self.config = config

        # select general or electrode dict
        if trode is None:
            values = self.values
            func_args = ()
        else:
            values = self.values[trode]
            # the electrode should be given to the function as argument
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
        """Ambipolar diffusivity
        """
        zp = self.config['zp']
        zm = self.config['zm']
        Dp = self.config['Dp']
        Dm = self.config['Dm']
        return ((zp - zm) * Dp * Dm) / (zp * Dp - zm * Dm)

    def tp(self):
        """Cation transference number
        """
        zp = self.config['zp']
        zm = self.config['zm']
        Dp = self.config['Dp']
        Dm = self.config['Dm']
        return zp * Dp / (zp * Dp - zm * Dm)

    def t_ref(self):
        """Reference time scale
        """
        return self.config['L_ref']**2 / self.config['D_ref']

    def curr_ref(self):
        """Reference current
        """
        return 3600. / self.config['t_ref']

    def sigma_s_ref(self):
        """Reference conductivity
        """
        return self.config['L_ref']**2 * constants.F**2 * constants.c_ref \
            / (self.config['t_ref'] * constants.k * constants.N_A * constants.T_ref)

    def currset(self):
        """Total current
        """
        return self.config['1C_current_density'] * self.config['Crate']  # A/m^2

    def Rser_ref(self):
        """Reference series resistance
        """
        return constants.k * constants.T_ref \
            / (constants.e * self.config['curr_ref'] * self.config['1C_current_density'])

    def csmax(self, trode):
        """Maximum concentration for given electrode
        """
        return self.config[trode, 'rho_s'] / constants.N_A

    def cap(self, trode):
        """Theoretical capacity of given electrode
        """
        return constants.e * self.config['L'][trode] * (1 - self.config['poros'][trode]) \
            * self.config['P_L'][trode] * self.config[trode, 'rho_s']

    def numsegments(self):
        """Number of segments in voltage/current profile
        """
        return len(self.config['segments'])

    def L_ref(self):
        """Reference length scale
        """
        return self.config['L']['c']

    def D_ref(self):
        """Reference diffusivity
        """
        if self.config['elyteModelType'] == 'dilute':
            return self.config['Damb']
        else:
            return getattr(props_elyte, self.config['SMset'])()[-1]

    def z(self):
        """Electrode capacity ratio
        """
        if 'a' in self.config['trodes']:
            return self.config['c', 'cap'] / self.config['a', 'cap']
        else:
            # flat plate anode with assumed infinite supply of metal
            return 0.

    def limtrode(self):
        """Capacity-limiting electrode
        """
        if self.config['z'] < 1:
            return 'c'
        else:
            return 'a'

    def cs_ref(self, trode):
        """Reference concentration
        """
        if self.config[trode, 'type'] in constants.one_var_types:
            prefac = 1
        elif self.config[trode, 'type'] in constants.two_var_types:
            prefac = .5
        return prefac * self.config[trode, 'csmax']

    def muR_ref(self, trode):
        """Reference chemical potential of given electrode
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

    def phiRef(self, trode):
        """Reference electrostatic potential of given electrode
        """
        return -self.config[trode, 'muR_ref'][0]

    def power_ref(self):
        """Reference power of the system
        """
        return constants.k*constants.T_ref * \
            self.config['curr_ref']*self.config["1C_current_density"]/constants.e
