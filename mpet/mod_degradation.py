"""These models define degradation reactions in the particles,
only at the surface. Will be called from mod_electrodes selections
"""
import daetools.pyDAE as dae
import numpy as np

import mpet.ports as ports
import mpet.props_am as props_am
import mpet.electrode.reactions as reactions


class plating_none(dae.daeModel):
    def __init__(self, config, trode, vInd, pInd,
                 Name, Parent=None, Description=""):
        super().__init__(Name, Parent, Description)

        # Variables
        self.dcplatingbardt = dae.daeVariable("dcplatingbardt", dae.no_t, self,
                                              "Rate of plating growth on particle volume basis")
        # Ports
        self.portInPart = ports.portFromPart(
            "portInPart", dae.eInletPort, self,
            "Inlet port from electrolyte")
        self.c_lyte = self.portInPart.c_lyte
        self.phi_lyte = self.portInPart.phi_lyte
        self.phi_m = self.portInPart.phi_m

    def DeclareEquations(self):
        dae.daeModel.DeclareEquations(self)
        eq = self.CreateEquation("dcplatingbardt")
        eq.Residual = self.dcplatingbardt()  # no degradation

        for eq in self.Equations:
            eq.CheckUnitsConsistency = False


class plating_simple(dae.daeModel):
    """Gao et al., Joule 2020"""

    def __init__(self, config, trode, vInd, pInd,
                 Name, Parent=None, Description=""):
        super().__init__(Name, Parent, Description)

        self.config = config
        self.trode = trode
        self.ind = (vInd, pInd)

        # Variables
        self.dcplatingbardt = dae.daeVariable("dcplatingbardt", dae.no_t, self,
                                              "Rate of plating growth on particle volume basis")
        # Reaction rates of plating
        self.Rxn_pl = dae.daeVariable("Rxn_pl", dae.no_t, self, "Plating reaction rate")
        self.V_Li = dae.daeVariable("V_Li", dae.no_t, self, "Plated Li volume")

        # get SEI and plating reaction rates
        if config[trode, "Li_plating"]:
            self.calc_rxn_rate_plating = getattr(reactions,"plating_simple")

        # Ports
        self.portInPart = ports.portFromPart(
            "portInPart", dae.eInletPort, self,
            "Inlet port from electrolyte")

        self.c_lyte = self.portInPart.c_lyte
        self.phi_lyte = self.portInPart.phi_lyte
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

        muR_pl = calc_muR_plating(self.config, self.trode, self.ind, self.ISfuncs)
        eta_pl = calc_eta(muR_pl, mu_O)
        Rxn_pl = self.calc_rxn_rate_plating(self.V_Li(), eta_pl, self.get_trode_param("k0_pl"),
                                            T, self.get_trode_param("alpha_pl"))
        eq = self.CreateEquation("Rxn_plating")
        eq.Residual = self.Rxn_pl() - Rxn_pl[0]

        eq = self.CreateEquation("Lithium_plating_growth")
        eq.Residual = self.get_trode_param("psd_area") * self.Rxn_pl() - self.V_Li.dt() \
            / self.get_trode_param("Li_mm")

        eq = self.CreateEquation("Plating_current")
        eq.Residual = self.dcplatingbardt() - self.get_trode_param("delta_L")*self.Rxn_pl()

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


def calc_muR_plating(config, trode, ind, ISfuncs=None):
    muRfunc = props_am.muRfuncs(config, trode, ind).muR_pl
    muR_ref = config[trode, "muR_ref"]
    muR, actR = muRfunc(0, 0, muR_ref, ISfuncs)
    return muR
