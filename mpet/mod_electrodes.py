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
import scipy.special as spcl

import mpet.extern_funcs as extern_funcs
import mpet.geometry as geo
import mpet.ports as ports
import mpet.props_am as props_am
import mpet.props_elyte as props_elyte
import mpet.utils as utils
import mpet.electrode.reactions as reactions
from mpet.daeVariableTypes import *

class Mod2var(dae.daeModel):
    def __init__(self, Name, Parent=None, Description="", ndD=None,
                 ndD_s=None):
        dae.daeModel.__init__(self, Name, Parent, Description)
        if (ndD is None) or (ndD_s is None):
            raise Exception("Need input parameter dictionary")
        self.ndD = ndD
        self.ndD_s = ndD_s

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
        if ndD["type"] not in ["ACR2"]:
            self.Rxn1 = dae.daeVariable("Rxn1", dae.no_t, self, "Rate of reaction 1")
            self.Rxn2 = dae.daeVariable("Rxn2", dae.no_t, self, "Rate of reaction 2")
        else:
            self.Rxn1 = dae.daeVariable("Rxn1", dae.no_t, self, "Rate of reaction 1", [self.Dmn])
            self.Rxn2 = dae.daeVariable("Rxn2", dae.no_t, self, "Rate of reaction 2", [self.Dmn])

        #Get reaction rate function from dictionary name
        self.calc_rxn_rate=getattr(reactions,ndD["rxnType"])

        # Ports
        self.portInLyte = ports.portFromElyte(
            "portInLyte", dae.eInletPort, self, "Inlet port from electrolyte")
        self.portInBulk = ports.portFromBulk(
            "portInBulk", dae.eInletPort, self,
            "Inlet port from e- conducting phase")
        self.portInSystem = ports.portFromSystem(
            "portInSystem", dae.eInletPort, self,
            "Inlet port from charge/discharge")
        self.phi_lyte = self.portInLyte.phi_lyte
        self.c_lyte = self.portInLyte.c_lyte
        self.phi_m = self.portInBulk.phi_m
        self.charge_discharge = self.portInSystem.charge_discharge
        self.cycle_number = self.portInSystem.cycle_number

    def DeclareEquations(self):
        dae.daeModel.DeclareEquations(self)
        ndD = self.ndD
        N = ndD["N"]  # number of grid points in particle
        T = self.ndD_s["T"]  # nondimensional temperature
        r_vec, volfrac_vec = geo.get_unit_solid_discr(ndD['shape'], N)

        # Prepare the Ideal Solution log ratio terms
        self.ISfuncs1 = self.ISfuncs2 = None
        if ndD["logPad"]:
            self.ISfuncs1 = np.array([
                extern_funcs.LogRatio("LR1", self, dae.unit(), self.c1(k))
                for k in range(N)])
            self.ISfuncs2 = np.array([
                extern_funcs.LogRatio("LR2", self, dae.unit(), self.c2(k))
                for k in range(N)])
        ISfuncs = (self.ISfuncs1, self.ISfuncs2)

        # Prepare noise
        self.noise1 = self.noise2 = None
        if ndD["noise"]:
            numnoise = ndD["numnoise"]
            noise_prefac = ndD["noise_prefac"]
            tvec = np.linspace(0., 1.05*self.ndD_s["tend"], numnoise)
            noise_data1 = noise_prefac*np.random.randn(numnoise, N)
            noise_data2 = noise_prefac*np.random.randn(numnoise, N)
            # Previous_output is common for all external functions
            previous_output1 = []
            previous_output2 = []
            self.noise1 = [extern_funcs.InterpTimeVector(
                "noise1", self, dae.unit(), dae.Time(), tvec,
                noise_data1, previous_output1, _position_)
                for _position_ in range(N)]
            self.noise2 = [extern_funcs.InterpTimeVector(
                "noise2", self, dae.unit(), dae.Time(), tvec,
                noise_data2, previous_output2, _position_)
                for _position_ in range(N)]
        noises = (self.noise1, self.noise2)

        # Figure out mu_O, mu of the oxidized state
        mu_O, act_lyte = calc_mu_O(
            self.c_lyte(), self.phi_lyte(), self.phi_m(), T,
            self.ndD_s["elyteModelType"])

        # Define average filling fractions in particle
        eq1 = self.CreateEquation("c1bar")
        eq2 = self.CreateEquation("c2bar")
        eq1.Residual = self.c1bar()
        eq2.Residual = self.c2bar()
        for k in range(N):
            eq1.Residual -= self.c1(k) * volfrac_vec[k]
            eq2.Residual -= self.c2(k) * volfrac_vec[k]
        eq = self.CreateEquation("cbar")
        eq.Residual = self.cbar() - utils.mean_linear(self.c1bar(), self.c2bar())

        # Define average rate of filling of particle
        eq = self.CreateEquation("dcbardt")
        eq.Residual = self.dcbardt()
        for k in range(N):
            eq.Residual -= utils.mean_linear(self.c1.dt(k), self.c2.dt(k)) * volfrac_vec[k]

        c1 = np.empty(N, dtype=object)
        c2 = np.empty(N, dtype=object)
        c1[:] = [self.c1(k) for k in range(N)]
        c2[:] = [self.c2(k) for k in range(N)]
        if ndD["type"] in ["diffn2", "CHR2"]:
            # Equations for 1D particles of 1 field varible
            self.sld_dynamics_1D2var(c1, c2, mu_O, act_lyte, ISfuncs, noises)
        elif ndD["type"] in ["homog2", "homog2_sdn"]:
            # Equations for 0D particles of 1 field variables
            self.sld_dynamics_0D2var(c1, c2, mu_O, act_lyte, ISfuncs, noises)

        for eq in self.Equations:
            eq.CheckUnitsConsistency = False

    def sld_dynamics_0D2var(self, c1, c2, muO, act_lyte, ISfuncs, noises):
        ndD = self.ndD
        T = self.ndD_s["T"]
        c1_surf = c1
        c2_surf = c2
        (mu1R_surf, mu2R_surf), (act1R_surf, act2R_surf) = (calc_muR(
            (c1_surf, c2_surf), (self.c1bar(), self.c2bar()), T, ndD, ISfuncs))
        eta1 = calc_eta(mu1R_surf, muO)
        eta2 = calc_eta(mu2R_surf, muO)
        eta1_eff = eta1 + self.Rxn1()*ndD["Rfilm"]
        eta2_eff = eta2 + self.Rxn2()*ndD["Rfilm"]
        Rxn1 = self.calc_rxn_rate(
            eta1_eff, c1_surf, self.c_lyte(), ndD["k0"], T,
            act1R_surf, act_lyte, ndD["lambda"], ndD["alpha"])
        Rxn2 = self.calc_rxn_rate(
            eta2_eff, c2_surf, self.c_lyte(), ndD["k0"], T,
            act2R_surf, act_lyte, ndD["lambda"], ndD["alpha"])
        eq1 = self.CreateEquation("Rxn1")
        eq2 = self.CreateEquation("Rxn2")
        eq1.Residual = self.Rxn1() - Rxn1[0]
        eq2.Residual = self.Rxn2() - Rxn2[0]

        noise1, noise2 = noises
        eq1 = self.CreateEquation("dc1sdt")
        eq2 = self.CreateEquation("dc2sdt")
        eq1.Residual = self.c1.dt(0) - ndD["delta_L"]*Rxn1[0]
        eq2.Residual = self.c2.dt(0) - ndD["delta_L"]*Rxn2[0]
        if ndD["noise"]:
            eq1.Residual += noise1[0]()
            eq2.Residual += noise2[0]()

    def sld_dynamics_1D2var(self, c1, c2, muO, act_lyte, ISfuncs, noises):
        ndD = self.ndD
        N = ndD["N"]
        T = self.ndD_s["T"]
        # Equations for concentration evolution
        # Mass matrix, M, where M*dcdt = RHS, where c and RHS are vectors
        Mmat = get_Mmat(ndD['shape'], N)
        dr, edges = geo.get_dr_edges(ndD['shape'], N)

        # Get solid particle chemical potential, overpotential, reaction rate
        if ndD["type"] in ["diffn2", "CHR2"]:
            (mu1R, mu2R), (act1R, act2R) = calc_muR(
                (c1, c2), (self.c1bar(), self.c2bar()), T, ndD, ISfuncs)
            c1_surf = c1[-1]
            c2_surf = c2[-1]
            mu1R_surf, act1R_surf = mu1R[-1], act1R[-1]
            mu2R_surf, act2R_surf = mu2R[-1], act2R[-1]
        eta1 = calc_eta(mu1R_surf, muO)
        eta2 = calc_eta(mu2R_surf, muO)
        if ndD["type"] in ["ACR2"]:
            eta1_eff = np.array([eta1[i] + self.Rxn1(i)*ndD["Rfilm"] for i in range(N)])
            eta2_eff = np.array([eta2[i] + self.Rxn2(i)*ndD["Rfilm"] for i in range(N)])
        else:
            eta1_eff = eta1 + self.Rxn1()*ndD["Rfilm"]
            eta2_eff = eta2 + self.Rxn2()*ndD["Rfilm"]
        Rxn1 = self.calc_rxn_rate(
            eta1_eff, c1_surf, self.c_lyte(), ndD["k0"], T,
            act1R_surf, act_lyte, ndD["lambda"], ndD["alpha"])
        Rxn2 = self.calc_rxn_rate(
            eta2_eff, c2_surf, self.c_lyte(), ndD["k0"], T,
            act2R_surf, act_lyte, ndD["lambda"], ndD["alpha"])
        if ndD["type"] in ["ACR2"]:
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
        if ndD["type"] in ["ACR2"]:
            #since ACR2 has N equations, use the last N vars and drop the first 
            RHS1 = np.array([ndD["delta_L"]*self.Rxn1(i) for i in range(N)])
            RHS2 = np.array([ndD["delta_L"]*self.Rxn2(i) for i in range(N)])
        elif ndD["type"] in ["diffn2", "CHR2"]:
            # Positive reaction (reduction, intercalation) is negative
            # flux of Li at the surface.
            Flux1_bc = -0.5 * self.Rxn1()
            Flux2_bc = -0.5 * self.Rxn2()
            Dfunc = props_am.Dfuncs(ndD["Dfunc"]).Dfunc
            if ndD["type"] == "diffn2":
                pass
#                Flux1_vec, Flux2_vec = calc_Flux_diffn2(
#                    c1, c2, ndD["D"], Flux1_bc, Flux2_bc, dr, T)
            elif ndD["type"] == "CHR2":
                Flux1_vec, Flux2_vec = calc_flux_CHR2(
                    c1, c2, mu1R, mu2R, ndD["D"], Dfunc, Flux1_bc, Flux2_bc, dr, T)
            if ndD["shape"] == "sphere":
                area_vec = 4*np.pi*edges**2
            elif ndD["shape"] == "cylinder":
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
        noise1, noise2 = noises
        for k in range(N):
            eq1 = self.CreateEquation("dc1sdt_discr{k}".format(k=k))
            eq2 = self.CreateEquation("dc2sdt_discr{k}".format(k=k))
            eq1.Residual = LHS1_vec[k] - RHS1[k]
            eq2.Residual = LHS2_vec[k] - RHS2[k]
            if ndD["noise"]:
                eq1.Residual += noise1[k]()
                eq2.Residual += noise2[k]()


class Mod1var(dae.daeModel):
    def __init__(self, Name, Parent=None, Description="", ndD=None,
                 ndD_s=None):
        dae.daeModel.__init__(self, Name, Parent, Description)

        if (ndD is None) or (ndD_s is None):
            raise Exception("Need input parameter dictionary")
        self.ndD = ndD
        self.ndD_s = ndD_s

        # Domain
        self.Dmn = dae.daeDomain("discretizationDomain", self, dae.unit(),
                                 "discretization domain")

        # Variables
        self.c = dae.daeVariable("c", mole_frac_t, self,
                                 "Concentration in active particle",
                                 [self.Dmn])
        self.cbar = dae.daeVariable(
            "cbar", mole_frac_t, self,
            "Average concentration in active particle")
        self.dcbardt = dae.daeVariable("dcbardt", dae.no_t, self, "Rate of particle filling")
        self.dcSEIbardt = dae.daeVariable("dcSEIbardt", dae.no_t, self,
                "Rate of SEI growth on particle volume basis")
        if ndD["type"] not in ["ACR"]:
            self.Rxn = dae.daeVariable("Rxn", dae.no_t, self, "Rate of reaction")
            self.Rxn_sign = dae.daeVariable("Rxn_sign", dae.no_t, self, "Rate of reaction is +/-")
        else:
            self.Rxn = dae.daeVariable("Rxn", dae.no_t, self, "Rate of reaction", [self.Dmn])

        self.Rxn_SEI = dae.daeVariable("Rxn_SEI", dae.no_t, self,
                "Rate of SEI growth reaction")
        self.L1 = dae.daeVariable("L1", dae.no_t, self,
                "Primary SEI thickness")
        self.L2 = dae.daeVariable("L2", dae.no_t, self,
                "Secondary SEI thickness")
        self.c_eff_lyte = dae.daeVariable("c_eff_lyte", dae.no_t, self,
                "Lithium ion concentration in the primary SEI layer")
        self.c_solv = dae.daeVariable("c_solv", mole_frac_t, self, "Solvent concentration")
        self.phi_SEI_L1 = dae.daeVariable("phi_SEI_L1", elec_pot_t, self,
                "Electrostatic potential at SEI/electrolyte interface")
        self.phi_SEI_L0 = dae.daeVariable("phi_SEI_L0", elec_pot_t, self,
                "Electrostatic potential at electrode/SEI interface")

        #get SEI and plating reaction rates
        if ndD["SEI"]:
            self.calc_rxn_rate_SEI = getattr(reactions,"SEI")
        if ndD["Plating"]:
            self.calc_rxn_rate_plating = getattr(reactions,"plating_simple")

        #Get reaction rate function from dictionary name
        self.calc_rxn_rate=getattr(reactions,ndD["rxnType"])

        # Ports
        self.portInLyte = ports.portFromElyte(
            "portInLyte", dae.eInletPort, self,
            "Inlet port from electrolyte")
        self.portInBulk = ports.portFromBulk(
            "portInBulk", dae.eInletPort, self,
            "Inlet port from e- conducting phase")
        self.portInSystem = ports.portFromSystem(
            "portInSystem", dae.eInletPort, self,
            "Inlet port for charge_discharge") ### SD added 07/15/2020

        self.phi_lyte = self.portInLyte.phi_lyte
        self.c_lyte = self.portInLyte.c_lyte
        self.phi_m = self.portInBulk.phi_m
        self.charge_discharge = self.portInSystem.charge_discharge 
        self.cycle_number = self.portInSystem.cycle_number

    def DeclareEquations(self):
        dae.daeModel.DeclareEquations(self)
        ndD = self.ndD
        N = ndD["N"]  # number of grid points in particle
        T = self.ndD_s["T"]  # nondimensional temperature
        r_vec, volfrac_vec = geo.get_unit_solid_discr(ndD['shape'], N)

        # Prepare the Ideal Solution log ratio terms
        self.ISfuncs = None
        if ndD["logPad"]:
            self.ISfuncs = np.array([
                extern_funcs.LogRatio("LR", self, dae.unit(), self.c(k))
                for k in range(N)])

        # Prepare noise
        self.noise = None
        if ndD["noise"]:
            numnoise = ndD["numnoise"]
            noise_prefac = ndD["noise_prefac"]
            tvec = np.linspace(0., 1.05*self.ndD_s["tend"], numnoise)
            noise_data = noise_prefac*np.random.randn(numnoise, N)
            # Previous_output is common for all external functions
            previous_output = []
            self.noise = [
                extern_funcs.InterpTimeVector(
                    "noise", self, dae.unit(), dae.Time(), tvec, noise_data,
                    previous_output, _position_)
                for _position_ in range(N)]

        # calculate effective concentration from the Stern layer
        if ndD["SEI"]:

           R_SEI = ndD["R0SEI"]*(self.L1()/ndD["L10"])/(self.c_eff_lyte()**ndD["nu"])

           #resistance of the SEI layer to electrons
           eq = self.CreateEquation("resistance_SEI")
           eq.Residual = self.phi_SEI_L0() - self.phi_SEI_L1() - self.Rxn_SEI()*R_SEI

           #stern layer-double layer equilibrium between c_lyte and c_L1
           c_avg_1 = (self.c_eff_lyte()+self.c_lyte())/2
           #if concentrated solution model
           if self.ndD_s["elyteModelType"] == "SM":
               eq = self.CreateEquation("Stern_layer")
               thermFac, tp0 = getattr(props_elyte,self.ndD_s["SMset"])()[2:4]
               eq.Residual = self.c_eff_lyte() - self.c_lyte()*(np.exp(-(self.phi_lyte() - self.phi_SEI_L1())-ndD["E_ads"]))**(1/thermFac(c_avg_1))
           else:
               eq = self.CreateEquation("Stern_layer")
               eq.Residual = self.c_eff_lyte() - np.exp(-(self.phi_lyte() - self.phi_SEI_L1())-ndD["E_ads"])*self.c_lyte()

           #ionic flux between lithium 
           eq = self.CreateEquation("lithium_ion_reaction")
           eq.Residual = self.Rxn()*ndD["R0SEILi"] - (self.phi_SEI_L1()-self.phi_SEI_L0())

        else:

           eq = self.CreateEquation("resistance_SEI")
           eq.Residual = self.phi_SEI_L0() - self.phi_SEI_L1()

           eq = self.CreateEquation("resistance_SEI_1")
           eq.Residual = self.phi_lyte() - self.phi_SEI_L1()

           eq = self.CreateEquation("lithium_ion_reaction")
           eq.Residual = self.c_eff_lyte() - self.c_lyte()
           

        # F#igure out mu_O, mu of the oxidized state
        mu_O, act_lyte = calc_mu_O(self.c_eff_lyte(), self.phi_SEI_L0(), self.phi_m(), T,
                                   self.ndD_s["elyteModelType"])
        mu_O_SEI, act_lyte_SEI = calc_mu_O(self.c_eff_lyte(), self.phi_SEI_L1(), self.phi_m(), T,
                                   self.ndD_s["elyteModelType"])

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
        c[:] = [self.c(k) for k in range(N)]
        if ndD["type"] in ["ACR", "diffn", "CHR"]:
            # Equations for 1D particles of 1 field varible
            self.sld_dynamics_1D1var(c, mu_O, act_lyte, self.ISfuncs, self.noise)
        elif ndD["type"] in ["homog", "homog_sdn"]:
            # Equations for 0D particles of 1 field variables
            self.sld_dynamics_0D1var(c, mu_O, mu_O_SEI, act_lyte, self.ISfuncs, self.noise)
         
        p1 = np.exp(-self.L1()/ndD["zeta"]) 
        p2 = 1-self.c_solv()/self.ndD_s["c0_solv"]
        w1 = p1/(p1 + p2) # normalized to conserve mass
        w2 = p2/(p1 + p2) # normalized to conserve mass

        eq = self.CreateEquation("Solvent_diffusion_Fick")
        eq.Residual = (self.ndD_s["c0_solv"] - self.c_solv())*(self.ndD_s["D_solv"]/self.L2()) - self.Rxn_SEI()

        eq = self.CreateEquation("Primary_SEI_growth")
        eq.Residual = w1*self.Rxn_SEI() - ndD["vfrac_1"]*self.L1.dt()

        eq = self.CreateEquation("Secondary_SEI_growth")
        eq.Residual = w2*self.Rxn_SEI() - ndD["vfrac_2"]*self.L2.dt()

        for eq in self.Equations:
            eq.CheckUnitsConsistency = False

    def sld_dynamics_0D1var(self, c, muO, muO_SEI, act_lyte, ISfuncs, noise):
        ndD = self.ndD
        T = self.ndD_s["T"]
        c_surf = c
        muR_surf, actR_surf = calc_muR(c_surf, self.cbar(), T, ndD, ISfuncs)
        eta = calc_eta(muR_surf, muO)
        eta_eff = eta + self.Rxn()*ndD["Rfilm"]
        Rxn = self.calc_rxn_rate(
            eta_eff, c_surf, self.c_eff_lyte(), ndD["k0"], T,
            actR_surf, act_lyte, ndD["lambda"], ndD["alpha"])
        eq = self.CreateEquation("Rxn")
        eq.Residual = self.Rxn() - Rxn[0]

        # add SEI equations
        if ndD["SEI"]:

            muR_SEI, actR_SEI = calc_muR_SEI(c_surf, self.cbar(), T, ndD, ISfuncs)
            eta_SEI = calc_eta(muR_SEI, muO_SEI)
            Rxn_SEI = self.calc_rxn_rate_SEI(eta_SEI, self.c_eff_lyte(), self.c_eff_lyte(), self.c_solv(), ndD["k0_SEI"], T, ndD["alpha_SEI"])
            eq = self.CreateEquation("Rxn_SEI")
            eq.Residual = self.Rxn_SEI() - Rxn_SEI #convert to Rxn_deg[0] if space dependent

        else:
            eq = self.CreateEquation("Rxn_SEI")
            eq.Residual = self.Rxn_SEI() - 0 #convert to Rxn_deg[0] if space dependent

        eq = self.CreateEquation("dcsSEIdt")
        eq.Residual = self.dcSEIbardt() - ndD["delta_L"]*self.Rxn_SEI()

        #the sign function doesn't work sadly...
        self.IF(self.Rxn() >= 0)
        eq = self.CreateEquation("sign_Rxn")
        eq.Residual = self.Rxn_sign() - 1
        self.ELSE()
        eq = self.CreateEquation("sign_Rxn")
        eq.Residual = self.Rxn_sign() + 1
        self.END_IF()
 
        eq = self.CreateEquation("dcsdt")
        eq.Residual = self.c.dt(0) - ndD["delta_L"]*(self.Rxn()-self.Rxn_sign()*self.Rxn_SEI())
        if ndD["noise"]:
            eq.Residual += noise[0]()

    def sld_dynamics_1D1var(self, c, muO, act_lyte, ISfuncs, noise):
        ndD = self.ndD
        N = ndD["N"]
        T = self.ndD_s["T"]
        # Equations for concentration evolution
        # Mass matrix, M, where M*dcdt = RHS, where c and RHS are vectors
        Mmat = get_Mmat(ndD['shape'], N)
        dr, edges = geo.get_dr_edges(ndD['shape'], N)

        # Get solid particle chemical potential, overpotential, reaction rate
        if ndD["type"] in ["ACR"]:
            c_surf = c
            muR_surf, actR_surf = calc_muR(
                c_surf, self.cbar(), T, ndD, ISfuncs)
        elif ndD["type"] in ["diffn", "CHR"]:
            muR, actR = calc_muR(c, self.cbar(), T, ndD, ISfuncs)
            c_surf = c[-1]
            muR_surf = muR[-1]
            if actR is None:
                actR_surf = None
            else:
                actR_surf = actR[-1]
        eta = calc_eta(muR_surf, muO)
        if ndD["type"] in ["ACR"]:
            eta_eff = np.array([eta[i] + self.Rxn(i)*ndD["Rfilm"] for i in range(N)])
        else:
            eta_eff = eta + self.Rxn()*ndD["Rfilm"]
        Rxn = self.calc_rxn_rate(
            eta_eff, c_surf, self.c_lyte(), ndD["k0"], T,
            actR_surf, act_lyte, ndD["lambda"], ndD["alpha"])
        if ndD["type"] in ["ACR"]:
            for i in range(N):
                eq = self.CreateEquation("Rxn_{i}".format(i=i))
                eq.Residual = self.Rxn(i) - Rxn[i]
        else:
            eq = self.CreateEquation("Rxn")
            eq.Residual = self.Rxn() - Rxn

        # Get solid particle fluxes (if any) and RHS
        if ndD["type"] in ["ACR"]:
            #since ACR has N equations, use the last N variables
            RHS = np.array([ndD["delta_L"]*self.Rxn(i) for i in range(N)])
        elif ndD["type"] in ["diffn", "CHR"]:
            # Positive reaction (reduction, intercalation) is negative
            # flux of Li at the surface.
            Flux_bc = -self.Rxn()
            Dfunc = props_am.Dfuncs(ndD["Dfunc"]).Dfunc
            if ndD["type"] == "diffn":
                Flux_vec = calc_flux_diffn(c, ndD["D"], Dfunc, Flux_bc, dr, T)
            elif ndD["type"] == "CHR":
                Flux_vec = calc_flux_CHR(c, muR, ndD["D"], Dfunc, Flux_bc, dr, T)
            if ndD["shape"] == "sphere":
                area_vec = 4*np.pi*edges**2
            elif ndD["shape"] == "cylinder":
                area_vec = 2*np.pi*edges  # per unit height
            #find the rate
            RHS = -np.diff(Flux_vec * area_vec)

        dcdt_vec = np.empty(N, dtype=object)
        dcdt_vec[0:N] = [self.c.dt(k) for k in range(N)]
        LHS_vec = MX(Mmat, dcdt_vec)
        for k in range(N):
            eq = self.CreateEquation("dcsdt_discr{k}".format(k=k))
            eq.Residual = LHS_vec[k] - RHS[k]
            if ndD["noise"]:
                eq.Residual += noise[k]()


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


def calc_flux_diffn(c, D, Dfunc, Flux_bc, dr, T):
    N = len(c)
    Flux_vec = np.empty(N+1, dtype=object)
    Flux_vec[0] = 0  # Symmetry at r=0
    Flux_vec[-1] = Flux_bc
    c_edges = utils.mean_harmonic(c)
    Flux_vec[1:N] = -D/T * Dfunc(c_edges) * np.diff(c)/dr
    return Flux_vec


def calc_flux_CHR(c, mu, D, Dfunc, Flux_bc, dr, T):
    N = len(c)
    Flux_vec = np.empty(N+1, dtype=object)
    Flux_vec[0] = 0  # Symmetry at r=0
    Flux_vec[-1] = Flux_bc
    c_edges = utils.mean_harmonic(c)
    Flux_vec[1:N] = -D/T * Dfunc(c_edges) * np.diff(mu)/dr
    return Flux_vec


def calc_flux_CHR2(c1, c2, mu1_R, mu2_R, D, Dfunc, Flux1_bc, Flux2_bc, dr, T):
    N = len(c1)
    Flux1_vec = np.empty(N+1, dtype=object)
    Flux2_vec = np.empty(N+1, dtype=object)
    Flux1_vec[0] = 0.  # symmetry at r=0
    Flux2_vec[0] = 0.  # symmetry at r=0
    Flux1_vec[-1] = Flux1_bc
    Flux2_vec[-1] = Flux2_bc
    c1_edges = utils.mean_harmonic(c1)
    c2_edges = utils.mean_harmonic(c2)
    Flux1_vec[1:N] = -D/T * Dfunc(c1_edges) * np.diff(mu1_R)/dr
    Flux2_vec[1:N] = -D/T * Dfunc(c2_edges) * np.diff(mu2_R)/dr
    return Flux1_vec, Flux2_vec


def calc_mu_O(c_lyte, phi_lyte, phi_sld, T, elyteModelType):
    if elyteModelType == "SM":
        mu_lyte = phi_lyte
        act_lyte = c_lyte
    elif elyteModelType == "dilute":
        act_lyte = c_lyte
        mu_lyte = T*np.log(act_lyte) + phi_lyte
    mu_O = mu_lyte - phi_sld
    return mu_O, act_lyte


def calc_muR(c, cbar, T, ndD, ISfuncs=None):
    muRfunc = props_am.muRfuncs(T, ndD).muRfunc
    muR_ref = ndD["muR_ref"]
    muR, actR = muRfunc(c, cbar, muR_ref, ISfuncs)
    return muR, actR


def calc_muR_SEI(c, cbar, T, ndD, ISfuncs=None):
    muRfunc = props_am.muRfuncs(T, ndD).muR_SEI
    muR_ref = ndD["muR_ref"][0]
    muR, actR = muRfunc(c, cbar, muR_ref, ISfuncs)
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

def lyte_fluxes(c_eff, phi_lyte2, phi_lyte1, dxd1, eps_o_tau, ndD):
    zp, zm, nup, num = ndD["zp"], ndD["zm"], ndD["nup"], ndD["num"]
    nu = nup + num
    T = ndD["T"]
    if ndD["elyteModelType"] == "dilute":
        Dp = eps_o_tau * ndD["Dp"]
        Dm = eps_o_tau * ndD["Dm"]
#        Np_edges_int = nup*(-Dp*np.diff(c_lyte)/dxd1
#                            - Dp*zp*c_edges_int*np.diff(phi_lyte)/dxd1)
        Nm_edges_int = num*(- Dm/T*zm*c_eff*(phi_lyte2-phi_lyte1)/dxd1)
       # i_edges_int = (-((nup*zp*Dp + num*zm*Dm)*(c_lyte2-c_lyte1)/dxd1)
       #                - (nup*zp**2*Dp + num*zm**2*Dm)/T*c_edges_int*(phi_lyte2-phi_lyte1)/dxd1)
#        i_edges_int = zp*Np_edges_int + zm*Nm_edges_int
    elif ndD["elyteModelType"] == "SM":
        D_fs, sigma_fs, thermFac, tp0 = getattr(props_elyte,ndD["SMset"])()[:-1]
        # modify the free solution transport properties for porous media

        def D(c):
            return eps_o_tau*D_fs(c)

        def sigma(c):
            return eps_o_tau*sigma_fs(c)
        sp, n = ndD["sp"], ndD["n_refTrode"]
        i_edges_int = -sigma(c_eff)/T * (
            (phi_lyte2-phi_lyte1)/dxd1
            )
        Nm_edges_int = num*((1./(num*zm)*(1-tp0(c_eff))*i_edges_int))
    return Nm_edges_int
