from string import Formatter as StringFormatter
from mpet import props_elyte

from mpet.exceptions import UnknownParameterError
from mpet.config import constants


class DerivedValues:
    def __init__(self):
        """
        """
        # to keep track of values that were already calculated
        # initialize with empty dicts for electrodes
        self.values = {'c': {}, 'a': {}}

        self.config = None

        self.equations = {'Damb': '(({zp} - {zm}) * {Dp} * {Dm}) / ({zp} * {Dp} - {zm} * {Dm})',
                          'tp': '{zp} * {Dp} / ({zp} * {Dp} - {zm} * {Dm})',
                          't_ref': '{L_ref}**2 / {D_ref}',
                          'curr_ref': '3600. / {t_ref}',
                          'sigma_s_ref': '{L_ref}**2 - constants.F**2 * constants.c_ref / '
                                '({t_ref} * constants.k * constants.N_A * constants.T_ref)'}

    def __repr__(self):
        return dict.__repr__(self.values)

    def get(self, config, item, trode=None):
        """
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
                values[item] = self._process_equation(self.equations[item])
            else:
                # get the method to calculate the value
                try:
                    func = getattr(self, item)
                except AttributeError:
                    raise UnknownParameterError(f"Unknown parameter: {item}")

                try:
                    values[item] = func(*args)
                except TypeError:
                    # TypeError occurs when calling arg-less func with args,
                    # or the other way around
                    raise Exception(f"Requested parameter {item} without electrode specification")

        return values[item]

    def _process_equation(self, equation):
        """
        """
        # ** operator does not call __getitem__ which would load the parameters,
        # so do this manually before calling **
        params = [item[1] for item in StringFormatter().parse(equation) if item[1] is not None]
        for param in params:
            self.config[param]
        # evaluate the equation
        return eval(equation.format(**self.config))

    def numsegments(self):
        """
        Number of segments
        """
        return len(self.config['segs'])

    def L_ref(self):
        """
        reference L
        TODO: with or without underscore?
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
        """
        if 'a' in self.config.trodes:
            return self.config['c', 'cap'] / self.config['a', 'cap']
        else:
            # flat plate anode with assumed infinite supply of metal
            return 0.

    def csmax(self, trode):
        """
        Maximum concentraion in electrode solids, mol/m^3
        """
        return self.config[trode, 'rho_s'] / constants.N_A

    def cs_ref(self, trode):
        """
        TODO: this parameter is just a scaling of 1 or .5 of cs_ref
        Perhaps just storing that scaling is enough?
        """
        if self.config[trode, 'type'] in constants.one_var_types:
            prefac = 1
        elif self.config[trode, 'type'] in constants.two_var_types:
            prefac = .5
        return prefac * self.config[trode, 'csmax']

    def cap(self, trode):
        """
        C / mˆ2
        """
        raise NotImplementedError('need acces to both sys and electrode config')
        return constants.e * self.config['L'][self.trode] \
            * (1 - self.config['poros'][self.trode]) \
            * self.config['P_L'][self.trode] * self.config['rho_s']
