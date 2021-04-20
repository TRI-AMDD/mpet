"""This defines the ports by which mod_cell interacts with mod_electrodes."""
import daetools.pyDAE as dae

from mpet.daeVariableTypes import *

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


class portFromSystem(dae.daePort):
    def __init__(self, Name, PortType, Model, Description=""):
        dae.daePort.__init__(self, Name, PortType, Model, Description)
        self.charge_discharge = dae.daeVariable(
            "charge_discharge", dae.no_t, self,
            "+1 indicates charge, and -1 indicates discharge")
        self.cycle_number = dae.daeVariable(
            "cycle_number", dae.no_t, self,
            "keeps track of which cycle number we are on in the mpet simulations")
