import daetools.pyDAE as dae
from mpet import ports
from mpet.daeVariableTypes import mole_frac_t, elec_pot_t


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

        # Domain
        self.Dmn = dae.daeDomain("discretizationDomain", self, dae.unit(),
                                 "discretization domain")

        # Variables
        self.c = dae.daeVariable("c", mole_frac_t, self,
                                 "Concentration in interface",
                                 [self.Dmn])

        self.phi = dae.daeVariable("phi", elec_pot_t, self,
                                   "Electrical potential in interface",
                                   [self.Dmn])

        # Ports
        self.portInLyte = ports.portFromElyte(
            "portInLyte", dae.eInletPort, self,
            "Inlet port from electrolyte")

        self.portOutLyte = ports.portFromElyte(
            "portOutLyte", dae.eOutletPort, self,
            "Electrolyte port from interface to particles")

    def DeclareEquations(self):
        super().DeclareEquations()
        Nvol = self.ndD_s["Nvol_i"]

        # Elyte to interface connection

        # set first grid point of interface equal to elyte
        eq = self.CreateEquation("c_lyte_to_interface")
        eq.Residual = self.portInLyte.c_lyte() - self.c(0)

        eq = self.CreateEquation("phi_lyte_to_interface")
        eq.Residual = self.portInLyte.phi_lyte() - self.phi(0)

        # Interface to particle connection through output port

        # last grid point of interface is output to particle
        eq = self.CreateEquation("c_interface_to_particle")
        eq.Residual = self.portOutLyte.c_lyte() - self.c(Nvol-1)

        eq = self.CreateEquation("phi_interface_to_particle")
        eq.Residual = self.portOutLyte.phi_lyte() - self.phi(Nvol-1)

        # Interface internals

        # fake interface region: set each grid point equal to previous grid point
        for n in range(1, Nvol):
            eq = self.CreateEquation(f"c_interface{n}")
            eq.Residual = self.c(n) - self.c(n-1)

            eq = self.CreateEquation(f"phi_interface{n}")
            eq.Residual = self.phi(n) - self.phi(n-1)
