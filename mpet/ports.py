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


class portFromSys(dae.daePort):
    def __init__(self, Name, PortType, Model, Description=""):
        dae.daePort.__init__(self, Name, PortType, Model, Description)
        self.current = dae.daeVariable(
            "current", elec_pot_t, self,
            "Current in system")
        self.endCondition = dae.daeVariable(
            "endCondition", elec_pot_t, self,
            "End condition in system")
        self.phi_applied = dae.daeVariable(
            "phi_applied", elec_pot_t, self,
            "Applied electric potential in system")
        self.ffrac_limtrode = dae.daeVariable(
            "ffrac_limtrode", elec_pot_t, self,
            "Filling fraction of limiting electrode in system")
