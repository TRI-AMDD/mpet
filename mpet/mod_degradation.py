"""These models define degradation reactions in the particles,
only at the surface. Will be called from mod_electrodes selections
"""
import daetools.pyDAE as dae
import numpy as np

import mpet.ports as ports
import mpet.props_am as props_am
import mpet.electrode.reactions as reactions
from mpet.daeVariableTypes import mole_frac_t


class SEI_none(dae.daeModel):
    def __init__(self, config, trode, vInd, pInd,
                 Name, Parent=None, Description=""):
        super().__init__(Name, Parent, Description)

        # Variables
        self.dcSEIbardt = dae.daeVariable("dcSEIbardt", dae.no_t, self,
                                          "Rate of SEI growth on particle volume basis")

        # ports
        self.portInPart = ports.portFromPart(
            "portInPart", dae.eInletPort, self,
            "Inlet port from electrolyte")
        self.phi_lyte = self.portInPart.phi_lyte
        self.c_lyte = self.portInPart.c_lyte
        self.phi_m = self.portInPart.phi_m

    def DeclareEquations(self):
        dae.daeModel.DeclareEquations(self)
        eq = self.CreateEquation("dcSEIbardt")
        eq.Residual = self.dcSEIbardt()  # no degradation

        for eq in self.Equations:
            eq.CheckUnitsConsistency = False


class SEI_adsorption(dae.daeModel):
    def __init__(self, config, trode, vInd, pInd,
                 Name, Parent=None, Description=""):
        super().__init__(Name, Parent, Description)

        self.config = config
        self.trode = trode
        self.ind = (vInd, pInd)

        # Variables
        self.dcSEIbardt = dae.daeVariable("dcSEIbardt", dae.no_t, self,
                                          "Rate of SEI growth on particle volume basis")
        # Get reaction rate function from dictionary name
        self.Rxn_SEI = dae.daeVariable("Rxn_SEI", dae.no_t, self,
                                       "Rate of SEI growth reaction")
        self.L1 = dae.daeVariable("L1", dae.no_t, self,
                                  "Primary SEI thickness")
        self.L2 = dae.daeVariable("L2", dae.no_t, self,
                                  "Secondary SEI thickness")
        self.c_solv = dae.daeVariable("c_solv", mole_frac_t, self, "Solvent concentration")
        self.a_e_SEI = dae.daeVariable("a_e_SEI", mole_frac_t, self,
                                       "Electron actitivty in SEI layer")

        # get SEI and plating reaction rates
        if config[trode, "SEI"]:
            self.calc_rxn_rate_SEI = getattr(reactions,"SEI")

        # ports
        self.portInPart = ports.portFromPart(
            "portInPart", dae.eInletPort, self,
            "Inlet port from electrolyte")
        self.phi_lyte = self.portInPart.phi_lyte
        self.c_lyte = self.portInPart.c_lyte
        self.phi_m = self.portInPart.phi_m

    def get_trode_param(self, item):
        """
        Shorthand to retrieve electrode-specific value
        """
        value = self.config[self.trode, item]
        # check if it is a particle-specific parameter
        if item in self.config.params_per_particle:
            value = value[self.ind]
        return value

    def DeclareEquations(self):
        dae.daeModel.DeclareEquations(self)
        T = self.config["T"]  # nondimensional temperature

        # Prepare the Ideal Solution log ratio terms
        self.ISfuncs = None

        # Figure out mu_O, mu of the oxidized state
        mu_O, act_lyte = calc_mu_O(self.c_lyte(), self.phi_lyte(), self.phi_m(), T,
                                   self.config["elyteModelType"])

        eq = self.CreateEquation("Solvent_diffusion_Fick")
        eq.Residual = (self.config["c0_solv"] - self.c_solv()) * \
            (self.config["D_solv"]/self.L2()) - self.Rxn_SEI()

        eq = self.CreateEquation("Primary_SEI_growth")
        eq.Residual = np.exp(-self.L1()/self.get_trode_param("zeta"))*self.Rxn_SEI() - \
            self.get_trode_param("c_SEI")*self.get_trode_param("vfrac_1")*self.L1.dt()

        eq = self.CreateEquation("Secondary_SEI_growth")
        eq.Residual = (1-np.exp(-self.L1()/self.get_trode_param("zeta")))*self.Rxn_SEI() - \
            self.get_trode_param("c_SEI")*self.get_trode_param("vfrac_2")*self.L2.dt()

        eq = self.CreateEquation("dcsSEIdt")
        eq.Residual = self.dcSEIbardt() - self.get_trode_param("delta_L")*self.Rxn_SEI()

        muR_SEI = calc_muR_SEI(self.config, self.trode, self.ind, self.ISfuncs)
        eta_SEI = calc_eta(muR_SEI, mu_O)
        Rxn_SEI = self.calc_rxn_rate_SEI(eta_SEI, self.a_e_SEI(), self.c_lyte(),
                                         self.c_solv(), self.get_trode_param("k0_SEI"),
                                         T, self.get_trode_param("alpha_SEI"))
        eq = self.CreateEquation("Rxn_SEI")
        eq.Residual = self.Rxn_SEI() - Rxn_SEI[0]  # convert to rxn_deg[0] if space dependent
        eq = self.CreateEquation("a_e_SEI")
        eq.Residual = self.a_e_SEI() - 1/(1+np.exp(-((self.phi_lyte()-self.phi_m())
                                                     - self.get_trode_param("eta_p"))))

        for eq in self.Equations:
            eq.CheckUnitsConsistency = False


def calc_eta(muR, muO):
    return muR - muO


def calc_mu_O(c_lyte, phi_lyte, phi_sld, T, elyteModelType):
    if elyteModelType == "SM":
        mu_lyte = phi_lyte
        act_lyte = c_lyte
    elif elyteModelType == "dilute":
        act_lyte = c_lyte
        mu_lyte = T*np.log(act_lyte) + phi_lyte
    mu_O = mu_lyte - phi_sld
    return mu_O, act_lyte


def calc_muR_SEI(config, trode, ind, ISfuncs=None):
    muRfunc = props_am.muRfuncs(config, trode, ind).muR_SEI
    muR_ref = config[trode, "muR_ref"]
    muR, actR = muRfunc(0, 0, muR_ref, ISfuncs)
    return muR
