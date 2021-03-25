from mpet import props_elyte

from mpet.config import constants


class Definitions:
    def __init__(self):
        """
        """
        # to keep track of values that were already calculated
        self.values = {}

    def get(self, item, params):
        """
        """
        if item not in self.values:
            # calculate the value
            try:
                func = getattr(self, item)
            except AttributeError:
                # No equation for the requested item
                return None

            self.values[item] = func(params)

        return self.values[item]


class DefinitionsSystem(Definitions):
    """
    Derived parameters for system
    """
    def numsegments(self, p):
        """
        Number of segments
        """
        return len(['segs'])

    def L_ref(self, p):
        """
        reference L
        TODO: with or without underscore?
        """
        return p['L']['c']

    def D_ref(self, p):
        """
        reference D
        """
        if p['elyteModelType'] == 'dilute':
            return p['Damb']
        else:
            return getattr(props_elyte, p['SMset'])()[-1]

    def Damb(self, p):
        """
        Ambipolar diffusivity
        """
        return ((p['zp'] - p['zm']) * p['Dp'] * p['Dm']) \
            / (p['zp'] * p['Dp'] - p['zm'] * p['Dm'])

    def tp(self, p):
        """
        Cation transference number
        """
        return p['zp'] * p['Dp'] / (p['zp'] * p['Dp'] - p['zm'] * p['Dm'])

    def t_ref(self, p):
        """
        Diffusive time scale
        """
        return p['L_ref']**2 / p['D_ref']

    def curr_ref(self, p):
        """
        """
        return 3600. / p['t_ref']

    def sigma_s_ref(self, p):
        """
        """
        return p['L_ref']**2 * constants.F**2 * constants.c_ref / p['L_ref']

    def z(self, p):
        """
        """
        raise NotImplementedError('need acces to both sys and electrode config')


class DefinitionsElectrode(Definitions):
    """
    Derived parameters for electrode
    TODO: some params need access to system and electrode config
    """
    def __init__(self, electrode):
        super().__init__()
        assert electrode in ['cathode', 'anode'], f"Invalid electrode: {electrode}"
        # internally, only the first character is used
        self.trode = electrode[0]

    def csmax(self, p):
        """
        Maximum concentraion in electrode solids, mol/m^3
        """
        return p['rho_s'] / constants.N_A

    def cs_ref(self, p):
        """
        TODO: this parameter is just a scaling of 1 or .5 of cs_ref
        Perhaps just storing that scaling is enough?
        """
        if p['type'] in constants.one_var_types:
            prefac = 1
        elif p['type'] in constants.two_var_types:
            prefac = .5
        return prefac * p['csmax']

    def cap(self, p):
        """
        C / mˆ2
        """
        raise NotImplementedError('need acces to both sys and electrode config')
        # return constants.e * p_system['L'][self.trode] * (1 - p_system['poros'][self.trode]) \
        #     * p_system['P_L'][self.trode] * p['rho_s']
