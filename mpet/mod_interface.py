import daetools.pyDAE as dae
from mpet import ports
# from mpet.daeVariableTypes import mole_frac_t, elec_pot_t


"""
Model for the interface between the electrolyte and active particles
"""


class InterfaceRegion(dae.daeModel):
    def __init__(self, Name, Parent=None, Description="", ndD=None,
                 ndD_s=None):
        super().__init__(Name, Parent, Description)
        if (ndD is None) or (ndD_s is None):
            raise Exception("Need input parameter dictionary")
        self.ndD = ndD
        self.ndD_s = ndD_s

        # Ports
        self.portInLyte = ports.portFromElyte(
            "portInLyte", dae.eInletPort, self,
            "Inlet port from electrolyte")

        self.portOutLyte = ports.portFromElyte(
            "portOutLyte", dae.eOutletPort, self,
            "Electrolyte port from interface to particles")
        # self.portInBulk = ports.portFromBulk(
        #     "portInBulk", dae.eInletPort, self,
        #     "Inlet port from e- conducting phase")
        # self.c_lyte = dae.daeVariable(
        #     "c_lyte", mole_frac_t, self,
        #     "Concentration in the electrolyte")
        # self.phi_lyte = dae.daeVariable(
        #     "phi_lyte", elec_pot_t, self,
        #     "Electric potential in the electrolyte")
        # self.phi_m = self.portInBulk.phi_m

    def DeclareEquations(self):
        super().DeclareEquations()

        # temporary interface region model: connect input to output
        # so interface does effectively nothing
        eq = self.CreateEquation("c_lyte_interface")
        eq.Residual = self.portInLyte.c_lyte() - self.portOutLyte.c_lyte()

        eq = self.CreateEquation("phi_lyte_interface")
        eq.Residual = self.portInLyte.phi_lyte() - self.portOutLyte.phi_lyte()
