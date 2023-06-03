"""These models define degradation reactions in the particles,
only at the surface. Will be called from mod_electrodes selections
"""
import daetools.pyDAE as dae
import numpy as np

import mpet.ports as ports
import mpet.props_am as props_am
import mpet.electrode.reactions as reactions
import mpet.utils as utils


class no_degradation(dae.daeModel):
    def __init__(self, config, trode, vInd, pInd,
                 Name, Parent=None, Description=""):
        super().__init__(Name, Parent, Description)

        # Variables
        self.c_tilde = dae.daeVariable("c_tilde", dae.no_t, self,
                                      "Rate of plating growth on particle volume basis")
        self.R_f = dae.daeVariable("R_f", dae.no_t, self,
                                   "Rate of plating growth on particle volume basis")
        self.dcdegradationdt = dae.daeVariable("dcdegradationdt", dae.no_t, self,
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
        eq = self.CreateEquation("ctilde")
        eq.Residual = self.c_tilde.dt()  # no degradation

        eq = self.CreateEquation("R_f")
        eq.Residual = self.R_f.dt()  # no degradation

        eq = self.CreateEquation("dcdegradationdt")
        eq.Residual = self.dcdegradationdt()  # no degradation

        for eq in self.Equations:
            eq.CheckUnitsConsistency = False


class simple_degradation(dae.daeModel):
    def __init__(self, config, trode, vInd, pInd,
                 Name, Parent=None, Description=""):
        super().__init__(Name, Parent, Description)

        self.config = config
        self.trode = trode
        self.ind = (vInd, pInd)

        # Variables
        self.c_tilde = dae.daeVariable("c_tilde", dae.no_t, self,
                                      "Rate of plating growth on particle volume basis")
        self.Rxn_c_tilde = dae.daeVariable("Rxn_c_tilde", dae.no_t, self,
                                      "Rate of plating growth on particle volume basis") 
        self.R_f = dae.daeVariable("R_f", dae.no_t, self,
                                   "rate of plating growth on particle volume basis")
        self.Rxn_R_f = dae.daeVariable("Rxn_R_f", dae.no_t, self,
                                      "Rate of plating growth on particle volume basis") 
        self.dcdegradationdt = dae.daeVariable("dcdegradationdt", dae.no_t, self,
                                      "Rate of plating growth on particle volume basis")
 
        if config[trode, "degradation"]:
            rxnType = "reaction_Tafel"
            self.calc_rxn_rate_Tafel = utils.import_function(config[trode, "rxnType_filename"],
                                                             rxnType,
                                                             f"mpet.electrode.reactions.{rxnType}")

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

        muR_c_tilde = calc_muR_c_tilde(self.config, self.trode, self.ind, self.ISfuncs)
        muR_R_f = calc_muR_R_f(self.config, self.trode, self.ind, self.ISfuncs)
        eta_c_tilde = calc_eta(muR_c_tilde, mu_O)
        eta_R_f = calc_eta(muR_R_f, mu_O)
        Rxn_c_tilde = self.calc_rxn_rate_Tafel(eta_c_tilde, self.get_trode_param("k0_c_tilde"),
                                               T)
        Rxn_R_f = self.calc_rxn_rate_Tafel(eta_R_f, self.get_trode_param("k0_R_f"),
                                           T)
 
        eq = self.CreateEquation("Rxn_c_tilde")
        eq.Residual = self.Rxn_c_tilde() - Rxn_c_tilde[0]

        eq = self.CreateEquation("Rxn_R_f")
        eq.Residual = self.Rxn_R_f() - Rxn_R_f[0]

        eq = self.CreateEquation("c_tilde")
        eq.Residual = self.c_tilde.dt() - 1/self.get_trode_param("delta_L")*self.Rxn_c_tilde()

        eq = self.CreateEquation("R_f")
        eq.Residual = self.R_f.dt() - self.get_trode_param("beta_deg")*self.Rxn_R_f()

        eq = self.CreateEquation("dcdegradationdt")
        eq.Residual = self.dcdegradationdt() - self.get_trode_param("delta_L")*(self.Rxn_R_f() + self.Rxn_c_tilde()) # no degradation

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


def calc_muR_c_tilde(config, trode, ind, ISfuncs=None):
    muRfunc = props_am.muRfuncs(config, trode, ind).muR_func_c_tilde
    muR_ref = config[trode, "muR_ref"]
    muR, actR = muRfunc(0, 0, muR_ref, ISfuncs)
    return muR

def calc_muR_R_f(config, trode, ind, ISfuncs=None):
    muRfunc = props_am.muRfuncs(config, trode, ind).muR_func_R_f
    muR_ref = config[trode, "muR_ref"]
    muR, actR = muRfunc(0, 0, muR_ref, ISfuncs)
    return muR
