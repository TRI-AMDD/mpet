"""These models define individual particles of active material.

This includes the equations for both 1-parameter models and 2-parameters models defining
 - mass conservation (concentration evolution)
 - reaction rate at the surface of the particles
In each model class it has options for different types of particles:
 - homogeneous
 - Fick-like diffusion
 - Cahn-Hilliard (with reaction boundary condition)
 - Allen-Cahn (with reaction throughout the particle)

These models can be instantiated from the mod_cell module to simulate various types of active
materials within a battery electrode.
"""


import daetools.pyDAE as dae

import numpy as np
import scipy.sparse as sprs
import scipy.interpolate as sintrp

import mpet.geometry as geo
import mpet.ports as ports
import mpet.props_am as props_am
import mpet.utils as utils
from mpet.daeVariableTypes import mole_frac_t


class Mod2var(dae.daeModel):
    def __init__(self, config, trode, vInd, pInd,
                 Name, Parent=None, Description=""):
        super().__init__(Name, Parent, Description)

        self.config = config
        self.trode = trode
        self.ind = (vInd, pInd)

        # Domain
        self.Dmn = dae.daeDomain("discretizationDomain", self, dae.unit(),
                                 "discretization domain")

        # Variables
        self.c1 = dae.daeVariable(
            "c1", mole_frac_t, self,
            "Concentration in 'layer' 1 of active particle", [self.Dmn])
        self.c2 = dae.daeVariable(
            "c2", mole_frac_t, self,
            "Concentration in 'layer' 2 of active particle", [self.Dmn])
        self.cbar = dae.daeVariable(
            "cbar", mole_frac_t, self,
            "Average concentration in active particle")
        self.c1bar = dae.daeVariable(
            "c1bar", mole_frac_t, self,
            "Average concentration in 'layer' 1 of active particle")
        self.c2bar = dae.daeVariable(
            "c2bar", mole_frac_t, self,
            "Average concentration in 'layer' 2 of active particle")
        self.dcbardt = dae.daeVariable("dcbardt", dae.no_t, self, "Rate of particle filling")
        if self.get_trode_param("type") not in ["ACR2"]:
            self.Rxn1 = dae.daeVariable("Rxn1", dae.no_t, self, "Rate of reaction 1")
            self.Rxn2 = dae.daeVariable("Rxn2", dae.no_t, self, "Rate of reaction 2")
        else:
            self.Rxn1 = dae.daeVariable("Rxn1", dae.no_t, self, "Rate of reaction 1", [self.Dmn])
            self.Rxn2 = dae.daeVariable("Rxn2", dae.no_t, self, "Rate of reaction 2", [self.Dmn])

        # Get reaction rate function
        rxnType = config[trode, "rxnType"]
        self.calc_rxn_rate = utils.import_function(config[trode, "rxnType_filename"],
                                                   rxnType,
                                                   f"mpet.electrode.reactions.{rxnType}")

        # Ports
        self.portInLyte = ports.portFromElyte(
            "portInLyte", dae.eInletPort, self, "Inlet port from electrolyte")
        self.portInBulk = ports.portFromBulk(
            "portInBulk", dae.eInletPort, self,
            "Inlet port from e- conducting phase")
        self.phi_lyte = self.portInLyte.phi_lyte
        self.c_lyte = self.portInLyte.c_lyte
        self.phi_m = self.portInBulk.phi_m
        if self.config[f"simInterface_{trode}"]:
            self.portOutParticle = ports.portFromParticle(
                "portOutParticle", dae.eOutletPort, self,
                "Outlet port from particle")

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
        N = self.get_trode_param("N")  # number of grid points in particle
        T = self.config["T"]  # nondimensional temperature
        r_vec, volfrac_vec = geo.get_unit_solid_discr(self.get_trode_param('shape'), N)

        # Prepare noise
        self.noise1 = self.noise2 = None
        if self.get_trode_param("noise"):
            numnoise = self.get_trode_param("numnoise")
            noise_prefac = self.get_trode_param("noise_prefac")
            tvec = np.linspace(0., 1.05*self.config["tend"], numnoise)
            noise_data1 = noise_prefac*np.random.randn(numnoise, N)
            noise_data2 = noise_prefac*np.random.randn(numnoise, N)
            self.noise1 = sintrp.interp1d(tvec, noise_data1, axis=0,
                                          bounds_error=False, fill_value=0.)
            self.noise2 = sintrp.interp1d(tvec, noise_data2, axis=0,
                                          bounds_error=False, fill_value=0.)
        noises = (self.noise1, self.noise2)

        # Figure out mu_O, mu of the oxidized state
        mu_O, act_lyte = calc_mu_O(
            self.c_lyte(), self.phi_lyte(), self.phi_m(), T,
            self.config, self.trode)

        # Define average filling fractions in particle
        eq1 = self.CreateEquation("c1bar")
        eq2 = self.CreateEquation("c2bar")
        eq1.Residual = self.c1bar()
        eq2.Residual = self.c2bar()
        for k in range(N):
            eq1.Residual -= self.c1(k) * volfrac_vec[k]
            eq2.Residual -= self.c2(k) * volfrac_vec[k]
        eq = self.CreateEquation("cbar")
        eq.Residual = self.cbar() - (0.5*self.c1bar() + 0.5*self.c2bar())

        # Define average rate of filling of particle
        eq = self.CreateEquation("dcbardt")
        eq.Residual = self.dcbardt()
        for k in range(N):
            eq.Residual -= (0.5*self.c1.dt(k) + 0.5*self.c2.dt(k)) * volfrac_vec[k]

        c1 = np.empty(N, dtype=object)
        c2 = np.empty(N, dtype=object)
        c1[:] = [self.c1(k) for k in range(N)]
        c2[:] = [self.c2(k) for k in range(N)]
        if self.get_trode_param("type") in ["diffn2", "CHR2", "ACR2"]:
            # Equations for 1D particles of 2 field variables
            self.sld_dynamics_1D2var(c1, c2, mu_O, act_lyte, noises)
        elif self.get_trode_param("type") in ["homog2", "homog2_sdn"]:
            # Equations for 0D particles of 2 field variables
            self.sld_dynamics_0D2var(c1, c2, mu_O, act_lyte, noises)

        for eq in self.Equations:
            eq.CheckUnitsConsistency = False

        if self.config[f"simInterface_{self.trode}"]:
            # Output dcbardt to interface region
            eq = self.CreateEquation("particle_to_interface_dcbardt")
            eq.Residual = self.portOutParticle.dcbardt() - self.dcbardt()

    def sld_dynamics_0D2var(self, c1, c2, muO, act_lyte, noises):
        T = self.config["T"]
        c1_surf = c1
        c2_surf = c2
        (mu1R_surf, mu2R_surf), (act1R_surf, act2R_surf) = calc_muR(
            (c1_surf, c2_surf), (self.c1bar(), self.c2bar()), self.config,
            self.trode, self.ind)
        eta1 = calc_eta(mu1R_surf, muO)
        eta2 = calc_eta(mu2R_surf, muO)
        eta1_eff = eta1 + self.Rxn1()*self.get_trode_param("Rfilm")
        eta2_eff = eta2 + self.Rxn2()*self.get_trode_param("Rfilm")
        noise1, noise2 = noises
        if self.get_trode_param("noise"):
            eta1_eff += noise1(dae.Time().Value)
            eta2_eff += noise2(dae.Time().Value)
        Rxn1 = self.calc_rxn_rate(
            eta1_eff, c1_surf, self.c_lyte(), self.get_trode_param("k0"),
            self.get_trode_param("E_A"), T, act1R_surf, act_lyte,
            self.get_trode_param("lambda"), self.get_trode_param("alpha"))
        Rxn2 = self.calc_rxn_rate(
            eta2_eff, c2_surf, self.c_lyte(), self.get_trode_param("k0"),
            self.get_trode_param("E_A"), T, act2R_surf, act_lyte,
            self.get_trode_param("lambda"), self.get_trode_param("alpha"))
        eq1 = self.CreateEquation("Rxn1")
        eq2 = self.CreateEquation("Rxn2")
        eq1.Residual = self.Rxn1() - Rxn1[0]
        eq2.Residual = self.Rxn2() - Rxn2[0]

        eq1 = self.CreateEquation("dc1sdt")
        eq2 = self.CreateEquation("dc2sdt")
        eq1.Residual = self.c1.dt(0) - self.get_trode_param("delta_L")*Rxn1[0]
        eq2.Residual = self.c2.dt(0) - self.get_trode_param("delta_L")*Rxn2[0]

    def sld_dynamics_1D2var(self, c1, c2, muO, act_lyte, noises):
        N = self.get_trode_param("N")
        T = self.config["T"]
        # Equations for concentration evolution
        # Mass matrix, M, where M*dcdt = RHS, where c and RHS are vectors
        Mmat = get_Mmat(self.get_trode_param('shape'), N)
        dr, edges = geo.get_dr_edges(self.get_trode_param('shape'), N)

        # Get solid particle chemical potential, overpotential, reaction rate
        if self.get_trode_param("type") in ["diffn2", "CHR2"]:
            (mu1R, mu2R), (act1R, act2R) = calc_muR(
                (c1, c2), (self.c1bar(), self.c2bar()), self.config, self.trode, self.ind)
            c1_surf = c1[-1]
            c2_surf = c2[-1]
            mu1R_surf, act1R_surf = mu1R[-1], act1R[-1]
            mu2R_surf, act2R_surf = mu2R[-1], act2R[-1]
            eta1 = calc_eta(mu1R_surf, muO)
            eta2 = calc_eta(mu2R_surf, muO)
        if self.get_trode_param("type") in ["ACR2"]:
            c1_surf = c1
            c2_surf = c2
            (mu1R, mu2R), (act1R, act2R) = calc_muR(
                (c1, c2), (self.c1bar(), self.c2bar()), self.config, self.trode, self.ind)
            mu1R_surf, act1R_surf = mu1R, act1R
            mu2R_surf, act2R_surf = mu2R, act2R
            eta1 = calc_eta(mu1R, muO)
            eta2 = calc_eta(mu2R, muO)
            eta1_eff = np.array([eta1[i]
                                 + self.Rxn1(i)*self.get_trode_param("Rfilm") for i in range(N)])
            eta2_eff = np.array([eta2[i]
                                 + self.Rxn2(i)*self.get_trode_param("Rfilm") for i in range(N)])
        else:
            eta1_eff = eta1 + self.Rxn1()*self.get_trode_param("Rfilm")
            eta2_eff = eta2 + self.Rxn2()*self.get_trode_param("Rfilm")
        Rxn1 = self.calc_rxn_rate(
            eta1_eff, c1_surf, self.c_lyte(), self.get_trode_param("k0"),
            self.get_trode_param("E_A"), T, act1R_surf, act_lyte,
            self.get_trode_param("lambda"), self.get_trode_param("alpha"))
        Rxn2 = self.calc_rxn_rate(
            eta2_eff, c2_surf, self.c_lyte(), self.get_trode_param("k0"),
            self.get_trode_param("E_A"), T, act2R_surf, act_lyte,
            self.get_trode_param("lambda"), self.get_trode_param("alpha"))
        if self.get_trode_param("type") in ["ACR2"]:
            for i in range(N):
                eq1 = self.CreateEquation("Rxn1_{i}".format(i=i))
                eq2 = self.CreateEquation("Rxn2_{i}".format(i=i))
                eq1.Residual = self.Rxn1(i) - Rxn1[i]
                eq2.Residual = self.Rxn2(i) - Rxn2[i]
        else:
            eq1 = self.CreateEquation("Rxn1")
            eq2 = self.CreateEquation("Rxn2")
            eq1.Residual = self.Rxn1() - Rxn1
            eq2.Residual = self.Rxn2() - Rxn2

        # Get solid particle fluxes (if any) and RHS
        if self.get_trode_param("type") in ["ACR2"]:
            RHS1 = 0.5*np.array([self.get_trode_param("delta_L")*self.Rxn1(i) for i in range(N)])
            RHS2 = 0.5*np.array([self.get_trode_param("delta_L")*self.Rxn2(i) for i in range(N)])
        elif self.get_trode_param("type") in ["diffn2", "CHR2"]:
            # Positive reaction (reduction, intercalation) is negative
            # flux of Li at the surface.
            Flux1_bc = -0.5 * self.Rxn1()
            Flux2_bc = -0.5 * self.Rxn2()
            Dfunc_name = self.get_trode_param("Dfunc")
            Dfunc = utils.import_function(self.get_trode_param("Dfunc_filename"),
                                          Dfunc_name,
                                          f"mpet.electrode.diffusion.{Dfunc_name}")
            if self.get_trode_param("type") == "CHR2":
                noise1, noise2 = noises
                Flux1_vec, Flux2_vec = calc_flux_CHR2(
                    c1, c2, mu1R, mu2R, self.get_trode_param("D"), Dfunc,
                    self.get_trode_param("E_D"), Flux1_bc, Flux2_bc, dr, T, noise1, noise2)
            if self.get_trode_param("shape") == "sphere":
                area_vec = 4*np.pi*edges**2
            elif self.get_trode_param("shape") == "cylinder":
                area_vec = 2*np.pi*edges  # per unit height
            RHS1 = -np.diff(Flux1_vec * area_vec)
            RHS2 = -np.diff(Flux2_vec * area_vec)
#            kinterlayer = 1e-3
#            interLayerRxn = (kinterlayer * (1 - c1) * (1 - c2) * (act1R - act2R))
#            RxnTerm1 = -interLayerRxn
#            RxnTerm2 = interLayerRxn
            RxnTerm1 = 0
            RxnTerm2 = 0
            RHS1 += RxnTerm1
            RHS2 += RxnTerm2

        dc1dt_vec = np.empty(N, dtype=object)
        dc2dt_vec = np.empty(N, dtype=object)
        dc1dt_vec[0:N] = [self.c1.dt(k) for k in range(N)]
        dc2dt_vec[0:N] = [self.c2.dt(k) for k in range(N)]
        LHS1_vec = MX(Mmat, dc1dt_vec)
        LHS2_vec = MX(Mmat, dc2dt_vec)
        for k in range(N):
            eq1 = self.CreateEquation("dc1sdt_discr{k}".format(k=k))
            eq2 = self.CreateEquation("dc2sdt_discr{k}".format(k=k))
            eq1.Residual = LHS1_vec[k] - RHS1[k]
            eq2.Residual = LHS2_vec[k] - RHS2[k]


class Mod1var(dae.daeModel):
    def __init__(self, config, trode, vInd, pInd,
                 Name, Parent=None, Description=""):
        super().__init__(Name, Parent, Description)

        self.config = config
        self.trode = trode
        self.ind = (vInd, pInd)

        # Domain
        self.Dmn = dae.daeDomain("discretizationDomain", self, dae.unit(),
                                 "discretization domain")
        # Variables
        self.c = dae.daeVariable("c", mole_frac_t, self,
                                 "Concentration in active particle",
                                 [self.Dmn])

        # Creation of the ghost points to assit BC

        if self.get_trode_param("type") in ["ACR_Diff"]:
            self.c_left_GP = dae.daeVariable("c_left", mole_frac_t, self,
                                             "Concentration on the left side of the particle")
            self.c_right_GP = dae.daeVariable("c_right", mole_frac_t, self,
                                              "Concentration on the right side of the particle")

        self.cbar = dae.daeVariable(
            "cbar", mole_frac_t, self,
            "Average concentration in active particle")
        self.dcbardt = dae.daeVariable("dcbardt", dae.no_t, self, "Rate of particle filling")
        if config[trode, "type"] not in ["ACR", "ACR_Diff"]:
            self.Rxn = dae.daeVariable("Rxn", dae.no_t, self, "Rate of reaction")
        else:
            self.Rxn = dae.daeVariable("Rxn", dae.no_t, self, "Rate of reaction", [self.Dmn])

        # Get reaction rate function
        rxnType = config[trode, "rxnType"]
        self.calc_rxn_rate = utils.import_function(config[trode, "rxnType_filename"],
                                                   rxnType,
                                                   f"mpet.electrode.reactions.{rxnType}")

        # Ports
        self.portInLyte = ports.portFromElyte(
            "portInLyte", dae.eInletPort, self,
            "Inlet port from electrolyte")
        self.portInBulk = ports.portFromBulk(
            "portInBulk", dae.eInletPort, self,
            "Inlet port from e- conducting phase")
        self.phi_lyte = self.portInLyte.phi_lyte
        self.c_lyte = self.portInLyte.c_lyte
        self.phi_m = self.portInBulk.phi_m
        if config[f"simInterface_{trode}"]:
            self.portOutParticle = ports.portFromParticle(
                "portOutParticle", dae.eOutletPort, self,
                "Outlet port from particle")

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
        N = self.get_trode_param("N")  # number of grid points in particle
        T = self.config["T"]  # nondimensional temperature
        r_vec, volfrac_vec = geo.get_unit_solid_discr(self.get_trode_param('shape'), N)

        # Prepare noise
        self.noise = None
        if self.get_trode_param("noise"):
            numnoise = self.get_trode_param("numnoise")
            noise_prefac = self.get_trode_param("noise_prefac")
            tvec = np.linspace(0., 1.05*self.config["tend"], numnoise)
            noise_data = noise_prefac*np.random.randn(numnoise, N)
            self.noise = sintrp.interp1d(tvec, noise_data, axis=0,
                                         bounds_error=False, fill_value=0.)

        # Figure out mu_O, mu of the oxidized state
        mu_O, act_lyte = calc_mu_O(self.c_lyte(), self.phi_lyte(), self.phi_m(), T,
                                   self.config, self.trode)

        # Define average filling fraction in particle
        eq = self.CreateEquation("cbar")
        eq.Residual = self.cbar()
        for k in range(N):
            eq.Residual -= self.c(k) * volfrac_vec[k]

        # Define average rate of filling of particle
        eq = self.CreateEquation("dcbardt")
        eq.Residual = self.dcbardt()
        for k in range(N):
            eq.Residual -= self.c.dt(k) * volfrac_vec[k]

        c = np.empty(N, dtype=object)
        if self.get_trode_param("type") in ["ACR_Diff"]:
            ctmp = utils.get_var_vec(self.c,N)
            c = np.hstack((self.c_left_GP(),ctmp,self.c_right_GP()))
        else:
            c[:] = [self.c(k) for k in range(N)]
        if self.get_trode_param("type") in ["ACR", "ACR_Diff", "diffn", "CHR"]:
            # Equations for 1D particles of 1 field varible
            self.sld_dynamics_1D1var(c, mu_O, act_lyte, self.noise)
        elif self.get_trode_param("type") in ["homog", "homog_sdn"]:
            # Equations for 0D particles of 1 field variables
            self.sld_dynamics_0D1var(c, mu_O, act_lyte, self.noise)

        for eq in self.Equations:
            eq.CheckUnitsConsistency = False

        if self.config[f"simInterface_{self.trode}"]:
            # Output dcbardt to interface region
            eq = self.CreateEquation("particle_to_interface_dcbardt")
            eq.Residual = self.portOutParticle.dcbardt() - self.dcbardt()

    def sld_dynamics_0D1var(self, c, muO, act_lyte, noise):
        T = self.config["T"]
        c_surf = c
        muR_surf, actR_surf = calc_muR(c_surf, self.cbar(), self.config,
                                       self.trode, self.ind)
        eta = calc_eta(muR_surf, muO)
        eta_eff = eta + self.Rxn()*self.get_trode_param("Rfilm")
        if self.get_trode_param("noise"):
            eta_eff += noise(dae.Time().Value)
        Rxn = self.calc_rxn_rate(
            eta_eff, c_surf, self.c_lyte(), self.get_trode_param("k0"),
            self.get_trode_param("E_A"), T, actR_surf, act_lyte,
            self.get_trode_param("lambda"), self.get_trode_param("alpha"))
        eq = self.CreateEquation("Rxn")
        eq.Residual = self.Rxn() - Rxn[0]

        eq = self.CreateEquation("dcsdt")
        eq.Residual = self.c.dt(0) - self.get_trode_param("delta_L")*self.Rxn()

    def sld_dynamics_1D1var(self, c, muO, act_lyte, noise):
        N = self.get_trode_param("N")
        T = self.config["T"]
        # Equations for concentration evolution
        # Mass matrix, M, where M*dcdt = RHS, where c and RHS are vectors
        Mmat = get_Mmat(self.get_trode_param('shape'), N)
        dr, edges = geo.get_dr_edges(self.get_trode_param('shape'), N)

        # Get solid particle chemical potential, overpotential, reaction rate
        if self.get_trode_param("type") in ["ACR", "ACR_Diff"]:
            c_surf = c
            # surface diffusion in the ACR C3 model
            if self.get_trode_param("type") in ["ACR_Diff"]:
                dx = 1/np.size(c)
                beta_s = self.get_trode_param("beta_s")
                eqL = self.CreateEquation("leftBC")
                eqL.Residual = (c_surf[0] - c_surf[1] +
                                - dx*beta_s*(c_surf[1]+0.008)*(1-c_surf[1]-0.008))
                eqR = self.CreateEquation("rightBC")
                eqR.Residual = (c_surf[-1] - c_surf[-2] +
                                - dx*beta_s*(c_surf[-2]+0.008)*(1-c_surf[-2]-0.008))

        if self.get_trode_param("type") in ["ACR", "ACR_Diff"]:
            muR_surf, actR_surf = calc_muR(
                c_surf, self.cbar(), self.config, self.trode, self.ind)
        elif self.get_trode_param("type") in ["diffn", "CHR"]:
            muR, actR = calc_muR(c, self.cbar(), self.config, self.trode, self.ind)
            c_surf = c[-1]
            muR_surf = muR[-1]
            if actR is None:
                actR_surf = None
            else:
                actR_surf = actR[-1]
        # surface diffusion in the ACR C3 model
        if self.get_trode_param("type") in ["ACR_Diff"]:
            eta = calc_eta(muR_surf[1:-1], muO)
        else:
            eta = calc_eta(muR_surf, muO)
        if self.get_trode_param("type") in ["ACR", "ACR_Diff"]:
            eta_eff = np.array([eta[i] + self.Rxn(i)*self.get_trode_param("Rfilm")
                                for i in range(N)])
        else:
            eta_eff = eta + self.Rxn()*self.get_trode_param("Rfilm")
        if self.get_trode_param("type") in ["ACR_Diff"]:
            Rxn = self.calc_rxn_rate(
                eta_eff, c_surf[1:-1], self.c_lyte(), self.get_trode_param("k0"),
                self.get_trode_param("E_A"), T, actR_surf[1:-1], act_lyte,
                self.get_trode_param("lambda"), self.get_trode_param("alpha"))
        else:
            Rxn = self.calc_rxn_rate(
                eta_eff, c_surf, self.c_lyte(), self.get_trode_param("k0"),
                self.get_trode_param("E_A"), T, actR_surf, act_lyte,
                self.get_trode_param("lambda"), self.get_trode_param("alpha"))
        if self.get_trode_param("type") in ["ACR", "ACR_Diff"]:
            for i in range(N):
                eq = self.CreateEquation("Rxn_{i}".format(i=i))
                eq.Residual = self.Rxn(i) - Rxn[i]
        else:
            eq = self.CreateEquation("Rxn")
            eq.Residual = self.Rxn() - Rxn

        # Get solid particle fluxes (if any) and RHS
        if self.get_trode_param("type") in ["ACR", "ACR_Diff"]:
            # For ACR model localized contact loss.
            if self.config['localized_losses']:
                gamma_con = self.get_trode_param('gamma_con')
                if gamma_con == 1:
                    RHS = np.array([self.get_trode_param("delta_L")*self.Rxn(i) for i in range(N)])
                else:
                    N_cont = int((gamma_con)*N)-12
                    RHS = np.zeros(N, dtype=object)
                    position = int(np.random.uniform()*N)
                    # random position of contact + 6 points on the sides to facilitate wetting
                    if 6+position+N_cont < N-6:
                        for i in range(0,6):
                            RHS[i] = self.get_trode_param("delta_L")*self.Rxn(i)
                        for i in range(N-6,N):
                            RHS[i] = self.get_trode_param("delta_L")*self.Rxn(i)
                        for i in range(position+6,position+N_cont,1):
                            RHS[i] = self.get_trode_param("delta_L")*self.Rxn(i)
                    else:
                        for i in range(0,6):
                            RHS[i] = self.get_trode_param("delta_L")*self.Rxn(i)
                        for i in range(position, N):
                            RHS[i] = self.get_trode_param("delta_L")*self.Rxn(i)
                        for i in range(6,N_cont-(N-position)):
                            RHS[i] = self.get_trode_param("delta_L")*self.Rxn(i)
            else:
                RHS = np.array([self.get_trode_param("delta_L")*self.Rxn(i) for i in range(N)])

        elif self.get_trode_param("type") in ["diffn", "CHR"]:
            # Positive reaction (reduction, intercalation) is negative
            # flux of Li at the surface.
            Flux_bc = -self.Rxn()
            Dfunc_name = self.get_trode_param("Dfunc")
            Dfunc = utils.import_function(self.get_trode_param("Dfunc_filename"),
                                          Dfunc_name,
                                          f"mpet.electrode.diffusion.{Dfunc_name}")
            if self.get_trode_param("type") == "diffn":
                Flux_vec = calc_flux_diffn(c, self.get_trode_param("D"), Dfunc,
                                           self.get_trode_param("E_D"), Flux_bc, dr, T, noise)
            elif self.get_trode_param("type") == "CHR":
                Flux_vec = calc_flux_CHR(c, muR, self.get_trode_param("D"), Dfunc,
                                         self.get_trode_param("E_D"), Flux_bc, dr, T, noise)
            if self.get_trode_param("shape") == "sphere":
                area_vec = 4*np.pi*edges**2
            elif self.get_trode_param("shape") == "cylinder":
                area_vec = 2*np.pi*edges  # per unit height
            RHS = -np.diff(Flux_vec * area_vec)

        dcdt_vec = np.empty(N, dtype=object)
        dcdt_vec[0:N] = [self.c.dt(k) for k in range(N)]
        LHS_vec = MX(Mmat, dcdt_vec)
        if self.get_trode_param("type") in ["ACR_Diff"]:
            # surface diffusion in the ACR C3 model
            surf_diff_vec = calc_surf_diff(c_surf, muR_surf, self.get_trode_param("D"))
        for k in range(N):
            eq = self.CreateEquation("dcsdt_discr{k}".format(k=k))
            if self.get_trode_param("type") in ["ACR_Diff"]:
                eq.Residual = LHS_vec[k] - RHS[k] - surf_diff_vec[k]
            else:
                eq.Residual = LHS_vec[k] - RHS[k]


# surface diffusion in the ACR C3 model
def calc_surf_diff(c_surf, muR_surf, D):
    N_2 = np.size(c_surf)
    dxs = 1./N_2
    c_surf_long = c_surf
    c_surf_short = (c_surf_long[0:-1] + c_surf_long[1:])/2
    surf_diff = D*(np.diff(c_surf_short*(1-c_surf_short)*np.diff(muR_surf)))/(dxs**2)
    return surf_diff


def calc_eta(muR, muO):
    return muR - muO


def get_Mmat(shape, N):
    r_vec, volfrac_vec = geo.get_unit_solid_discr(shape, N)
    if shape == "C3":
        Mmat = sprs.eye(N, N, format="csr")
    elif shape in ["sphere", "cylinder"]:
        Rs = 1.
        # For discretization background, see Zeng & Bazant 2013
        # Mass matrix is common for each shape, diffn or CHR
        if shape == "sphere":
            Vp = 4./3. * np.pi * Rs**3
        elif shape == "cylinder":
            Vp = np.pi * Rs**2  # per unit height
        vol_vec = Vp * volfrac_vec
        M1 = sprs.diags([1./8, 3./4, 1./8], [-1, 0, 1],
                        shape=(N, N), format="csr")
        M1[1,0] = M1[-2,-1] = 1./4
        M2 = sprs.diags(vol_vec, 0, format="csr")
        Mmat = M1*M2
    return Mmat


def calc_flux_diffn(c, D, Dfunc, E_D, Flux_bc, dr, T, noise):
    N = len(c)
    Flux_vec = np.empty(N+1, dtype=object)
    Flux_vec[0] = 0  # Symmetry at r=0
    Flux_vec[-1] = Flux_bc
    c_edges = utils.mean_linear(c)
    if noise is None:
        Flux_vec[1:N] = -D * Dfunc(c_edges) * np.exp(-E_D/T + E_D/1) * np.diff(c)/dr
    else:
        Flux_vec[1:N] = -D * Dfunc(c_edges) * np.exp(-E_D/T + E_D/1) * \
            np.diff(c + noise(dae.Time().Value))/dr
    return Flux_vec


def calc_flux_CHR(c, mu, D, Dfunc, E_D, Flux_bc, dr, T, noise):
    N = len(c)
    Flux_vec = np.empty(N+1, dtype=object)
    Flux_vec[0] = 0  # Symmetry at r=0
    Flux_vec[-1] = Flux_bc
    c_edges = utils.mean_linear(c)
    if noise is None:
        Flux_vec[1:N] = -D/T * Dfunc(c_edges) * np.exp(-E_D/T + E_D/1) * np.diff(mu)/dr
    else:
        Flux_vec[1:N] = -D/T * Dfunc(c_edges) * np.exp(-E_D/T + E_D/1) * \
            np.diff(mu + noise(dae.Time().Value))/dr
    return Flux_vec


def calc_flux_CHR2(c1, c2, mu1_R, mu2_R, D, Dfunc, E_D, Flux1_bc, Flux2_bc, dr, T, noise1, noise2):
    N = len(c1)
    Flux1_vec = np.empty(N+1, dtype=object)
    Flux2_vec = np.empty(N+1, dtype=object)
    Flux1_vec[0] = 0.  # symmetry at r=0
    Flux2_vec[0] = 0.  # symmetry at r=0
    Flux1_vec[-1] = Flux1_bc
    Flux2_vec[-1] = Flux2_bc
    c1_edges = utils.mean_linear(c1)
    c2_edges = utils.mean_linear(c2)
    if noise1 is None:
        Flux1_vec[1:N] = -D/T * Dfunc(c1_edges) * np.exp(-E_D/T + E_D/1) * np.diff(mu1_R)/dr
        Flux2_vec[1:N] = -D/T * Dfunc(c2_edges) * np.exp(-E_D/T + E_D/1) * np.diff(mu2_R)/dr
    else:
        Flux1_vec[1:N] = -D/T * Dfunc(c1_edges) * np.exp(-E_D/T + E_D/1) * \
            np.diff(mu1_R+noise1(dae.Time().Value))/dr
        Flux2_vec[1:N] = -D/T * Dfunc(c2_edges) * np.exp(-E_D/T + E_D/1) * \
            np.diff(mu2_R+noise2(dae.Time().Value))/dr
    return Flux1_vec, Flux2_vec


def calc_mu_O(c_lyte, phi_lyte, phi_sld, T, config, trode):
    elyteModelType = config["elyteModelType"]

    if config[f"simInterface_{trode}"]:
        elyteModelType = config["interfaceModelType"]

    if elyteModelType == "SM":
        mu_lyte = phi_lyte
        act_lyte = c_lyte
    elif elyteModelType == "dilute":
        act_lyte = c_lyte
        mu_lyte = T*np.log(act_lyte) + phi_lyte
    elif elyteModelType == "solid":
        a_slyte = config['a_slyte']
        cmax = config['cmax']
        act_lyte = (c_lyte / cmax) / (1 - c_lyte / cmax)*np.exp(a_slyte*(1 - 2*c_lyte))
        mu_lyte = phi_lyte
    mu_O = mu_lyte - phi_sld
    return mu_O, act_lyte


def calc_muR(c, cbar, config, trode, ind):
    muRfunc = props_am.muRfuncs(config, trode, ind).muRfunc
    muR_ref = config[trode, "muR_ref"]
    muR, actR = muRfunc(c, cbar, muR_ref)
    return muR, actR


def MX(mat, objvec):
    if not isinstance(mat, sprs.csr.csr_matrix):
        raise Exception("MX function designed for csr mult")
    n = objvec.shape[0]
    if isinstance(objvec[0], dae.pyCore.adouble):
        out = np.empty(n, dtype=object)
    else:
        out = np.zeros(n, dtype=float)
    # Loop through the rows
    for i in range(n):
        low = mat.indptr[i]
        up = mat.indptr[i+1]
        if up > low:
            out[i] = np.sum(
                mat.data[low:up] * objvec[mat.indices[low:up]])
        else:
            out[i] = 0.0
    return out
