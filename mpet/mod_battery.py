# import daetools.pyDAE as dae
# from pyUnits import s

# import numpy as np

# import mpet.extern_funcs as extern_funcs
# import mpet.geometry as geom
# import mpet.mod_electrodes as mod_electrodes
# import mpet.ports as ports
# import mpet.utils as utils
# from mpet.config import constants
# from mpet.daeVariableTypes import mole_frac_t, elec_pot_t, conc_t, temp_t

# # Dictionary of end conditions
# endConditions = {
#     1:"Vmax reached",
#     2:"Vmin reached"}

# class ModCell(dae.daeModel):
#     def __init__(self, config, Name, Parent=None, Description=""):
#         dae.daeModel.__init__(self, Name, Parent, Description)

#         self.config = config
#         Nlayers = config["Nlayers"]

#         self.DmnCell = dae.daeDomain(
#             "DmnCell", self, dae.unit(),
#             "Simulated layers in battery")
               
#         # Variables
#         self.phi_bound = {}
#         self.T_bound = {}

#         # ports
#         self.portsOutT = {}

#          for lInd in range(Nlayers):
#            self.portsOutLyte[lInd] = ports.portFromCell(
#                     "portCell{lInd}".format(lInd=lInd), dae.eOutletPort,
#                     self, "Electrolyte port to particles")