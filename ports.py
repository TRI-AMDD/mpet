"""This defines the ports by which mod_cell interacts with mod_electrodes."""
import daetools.pyDAE as dae

mole_frac_t = dae.daeVariableType(
    name="mole_frac_t", units=dae.unit(), lowerBound=0,
    upperBound=1, initialGuess=0.25, absTolerance=1e-6)
elec_pot_t = dae.daeVariableType(
    name="elec_pot_t", units=dae.unit(), lowerBound=-1e20,
    upperBound=1e20, initialGuess=0, absTolerance=1e-5)


class portFromElyte(dae.daePort):
    def __init__(self, Name, PortType, Model, Description=""):
        dae.daePort.__init__(self, Name, PortType, Model, Description)
        self.c_lyte = dae.daeVariable(
            "c_lyte", mole_frac_t, self,
            "Concentration in the electrolyte")
        self.phi_lyte = dae.daeVariable(
            "phi_lyte", elec_pot_t, self,
            "Electric potential in the electrolyte")


class portFromBulk(dae.daePort):
    def __init__(self, Name, PortType, Model, Description=""):
        dae.daePort.__init__(self, Name, PortType, Model, Description)
        self.phi_m = dae.daeVariable(
            "phi_m", elec_pot_t, self,
            "Electric potential in the e- conducting phase")
