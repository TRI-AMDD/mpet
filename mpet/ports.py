"""This defines the ports by which mod_cell interacts with mod_electrodes."""
import daetools.pyDAE as dae

from mpet.daeVariableTypes import mole_frac_t, elec_pot_t


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


class portFromParticle(dae.daePort):
    def __init__(self, Name, PortType, Model, Description=""):
        dae.daePort.__init__(self, Name, PortType, Model, Description)
        self.dcbardt = dae.daeVariable(
            "dcbardt", dae.no_t, self,
            "Rate of particle filling")


class portFromInterface(dae.daePort):
    def __init__(self, Name, PortType, Model, Description=""):
        dae.daePort.__init__(self, Name, PortType, Model, Description)
        self.i0 = dae.daeVariable(
            "i0", dae.no_t, self,
            "Current density in first interface region volume")

        self.Nm0 = dae.daeVariable(
            "Nm0", dae.no_t, self,
            "Flux in first interface region volume")
