from daetools.pyDAE import *

mole_frac_t = daeVariableType(name="mole_frac_t", units=unit(),
        lowerBound=0, upperBound=1, initialGuess=0.25,
        absTolerance=1e-6)
elec_pot_t = daeVariableType(name="elec_pot_t", units=unit(),
        lowerBound=-1e20, upperBound=1e20, initialGuess=0,
        absTolerance=1e-5)

class portFromElyte(daePort):
    def __init__(self, Name, PortType, Model, Description=""):
        daePort.__init__(self, Name, PortType, Model, Description)
        self.c_lyte = daeVariable("c_lyte", mole_frac_t, self,
                "Concentration in the electrolyte")
        self.phi_lyte = daeVariable("phi_lyte", elec_pot_t, self,
                "Electric potential in the electrolyte")

class portFromBulk(daePort):
    def __init__(self, Name, PortType, Model, Description=""):
        daePort.__init__(self, Name, PortType, Model, Description)
        self.phi_m = daeVariable("phi_m", elec_pot_t, self,
                "Electric potential in the e- conducting phase")
