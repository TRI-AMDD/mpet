import daetools.pyDAE as dae
# from pyUnits import s

# import numpy as np

# import mpet.extern_funcs as extern_funcs
# import mpet.geometry as geom
import mpet.mod_cell as mod_cell
import mpet.ports as ports
# import mpet.utils as utils
# from mpet.config import constants
from mpet.daeVariableTypes import mole_frac_t, elec_pot_t, temp_t

# Dictionary of end conditions
endConditions = {
    1:"Vmax reached",
    2:"Vmin reached"}


class ModBattery(dae.daeModel):
    def __init__(self, config, Name, Parent=None, Description=""):
        dae.daeModel.__init__(self, Name, Parent, Description)

        self.config = config
        Nlayers = config["Nlayers"]

        self.DmnBatt = dae.daeDomain(
            "DmnBatt", self, dae.unit(),
            "Simulated layers in battery")

        # Variables
        self.soc = dae.daeVariable(
            "soc", mole_frac_t, self,
            "Overall soc of a layer",
            [self.DmnBatt])
        self.T_battery = dae.daeVariable(
            "T_battery", temp_t, self,
            "Temperature along the battery",
            [self.DmnBatt])
        self.Q_layer = dae.daeVariable(
            "Q_layer", dae.no_t, self,
            "Rate of heat generation of each layer",
            [self.DmnBatt])
        self.I_layer = dae.daeVariable(
            "I_layer", dae.no_t, self,
            "Current passing through each layer",
            [self.DmnBatt])

        self.phi_applied = dae.daeVariable(
            "phi_applied", elec_pot_t, self,
            "Overall battery voltage")
        self.phi_cell = dae.daeVariable(
            "phi_cell", elec_pot_t, self,
            "Voltage between electrodes (phi_applied less series resistance)")
        self.current = dae.daeVariable(
            "current", dae.no_t, self, "Total current of the cell")
        self.endCondition = dae.daeVariable(
            "endCondition", dae.no_t, self, "A nonzero value halts the simulation")

        # ports
        self.portsOutT = {}
        self.layers = {}
        for lInd in range(Nlayers):
            self.portsOutLayer[lInd] = ports.portFromCell(
                "portCell{lInd}".format(lInd=lInd), dae.eOutletPort,
                self, "Electrolyte port to particles")
            cMod = mod_cell.ModCell
            self.layers[lInd] = cMod(config, lInd,
                                     Name="layer{lInd}".format(lInd=lInd),
                                     Parent=self)

            self.ConnectPorts(self.portsCap,
                              self.layers[lInd].portInCap)

    def DeclareEquations(self):
        dae.daeModel.DeclareEquations(self)

        # Some values of domain lengths
        config = self.config
        Nlayers = config["Nlayers"]
        dl = 1/Nlayers
        # sum of all I_layers = I 

