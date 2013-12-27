#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import sys
import os
import ConfigParser
from time import localtime, strftime

import numpy as np
import scipy.sparse as sprs
import scipy.interpolate as sint
import scipy.io as sio

from daetools.pyDAE import *
from daetools.pyDAE.data_reporters import *
from daetools.solvers.superlu import pySuperLU
#from daetools.solvers.superlu_mt import pySuperLU_MT
from daetools.solvers.trilinos import pyTrilinos
#from daetools.solvers.intel_pardiso import pyIntelPardiso
from pyUnits import s
#from pyUnits import s, kg, m, K, Pa, mol, J, W

# Define some variable types
mole_frac_t = daeVariableType(name="mole_frac_t", units=unit(),
        lowerBound=0, upperBound=1, initialGuess=0.25,
        absTolerance=1e-6)
elec_pot_t = daeVariableType(name="elec_pot_t", units=unit(),
        lowerBound=-1e20, upperBound=1e20, initialGuess=0,
        absTolerance=1e-5)

class modMPET(daeModel):
    def __init__(self, Name, Parent=None, Description="",
            Ntrode=None, numpart=None, P=None, simSurfCathCond=False,
            profileType="CC"):
        daeModel.__init__(self, Name, Parent, Description)

        if (Ntrode is None) or (numpart is None):
            raise Exception("Need particle size distr. as input")
        self.P = P
        self.profileType = profileType

        # Domains where variables are distributed
        if Ntrode > 1: # If we have a separator
            self.Nsep = daeDomain("Nsep", self, unit(),
                    "Number of control volumes in the separator")
        self.Ntrode = daeDomain("Ntrode", self, unit(),
                "Number of control volumes in the electrode")
        self.numpart = daeDomain("numpart", self, unit(),
                "Number of particles sampled per electrode control volume")
        self.Nsld_mat = np.empty((Ntrode, numpart), dtype=object)
        for i in range(Ntrode):
            for j in range(numpart):
                self.Nsld_mat[i, j] = daeDomain("vol{i}_part{j}".format(
                    i=i, j=j), self, unit(),
                    "Number of discretizations for particle "
                    + "j in volume i".format(i=i,j=j))

        # Variables
        self.c_lyte_trode = daeVariable("c_lyte_trode", mole_frac_t, self,
                "Concentration in the electrolyte in the electrode",
                [self.Ntrode])
        self.phi_lyte_trode = daeVariable("phi_lyte_trode", elec_pot_t, self,
                "Electrostatic potential in electrolyte in electrode",
                [self.Ntrode])
        if Ntrode > 1: # If we have a separator
            self.c_lyte_sep = daeVariable("c_lyte_sep", mole_frac_t, self,
                    "Concentration in the electrolyte in the separator",
                    [self.Nsep])
            self.phi_lyte_sep = daeVariable("phi_lyte_sep", elec_pot_t, self,
                    "Electrostatic potential in electrolyte in separator",
                    [self.Nsep])
        self.phi_applied = daeVariable("phi_applied", elec_pot_t, self,
                "Overall battery voltage (at anode current collector)")
        self.current = daeVariable("current", no_t, self,
                "Total current of the cell")
        self.c_sld = np.empty((Ntrode, numpart), dtype=object)
        for i in range(Ntrode):
            for j in range(numpart):
                self.c_sld[i, j] = daeVariable("solid_c_vol{i}_part{j}".format(
                    i=i, j=j), mole_frac_t, self,
                    "Concentration in each solid particle",
                    [self.Nsld_mat[i, j]])
        # Only make a variable of this if we have to -- it's a lot of
        # equations to keep track of for nothing if we don't need it.
        if simSurfCathCond:
            self.phi_sld = np.empty((Ntrode, numpart), dtype=object)
            for i in range(Ntrode):
                for j in range(numpart):
                    self.phi_sld[i, j] = daeVariable("solid_p_vol{i}_part{j}".format(
                        i=i, j=j), elec_pot_t, self,
                        "Electrostatic potential in each solid particle",
                        [self.Nsld_mat[i, j]])
        self.cbar_sld = daeVariable("cbar_sld", mole_frac_t, self,
                "Average concentration in each particle",
                [self.Ntrode, self.numpart])
        self.phi_c = daeVariable("phi_cath", elec_pot_t, self,
                "Electrostatic potential in the solid",
                [self.Ntrode])
        self.j_plus = daeVariable("j_plus", no_t, self,
                "Rate of reaction of positives per solid volume",
                [self.Ntrode])
        self.ffrac_cathode = daeVariable("ffrac_cathode",
                mole_frac_t, self,
                "Overall filling fraction of solids in cathode")

        # Parameters
        self.NumTrode = daeParameter("NumTrode", unit(), self,
                "Number of volumes in the electrode")
        self.NumPart = daeParameter("NumPart", unit(), self,
                "Number of particles in each electrode volume")
        self.NumSep = daeParameter("NumSep", unit(), self,
                "Number of volumes in the electrolyte")
        self.Lp = daeParameter("Lp", unit(), self,
                "loading percent (vol active per vol solid)")
        self.c_lyte0 = daeParameter("c_lyte0", unit(), self,
                "initial conc. of electrolyte [mol/m^3]")
        self.epsbeta = daeParameter("epsbeta", unit(), self,
                "porosity times beta")
        self.zp = daeParameter("zp", unit(), self,
                "cation charge number")
        self.zm = daeParameter("zm", unit(), self,
                "anion charge number")
        self.tp = daeParameter("tp", unit(), self,
                "positive transference number")
        self.poros_sep = daeParameter("poros_sep", unit(), self,
                "porosity in separator")
        self.poros_trode = daeParameter("poros_trode", unit(), self,
                "porosity in electrode")
        self.phi_cathode = daeParameter("phi_cathode", unit(), self,
                "potential at the cathode (phi_applied is relative to this)")
        self.td = daeParameter("td", unit(), self,
                "Diffusive time [s]")
        self.dim_Dp = daeParameter("dim_Dp", unit(), self,
                "diffusivity of positive ions [m^2/s]")
        self.dim_Dm = daeParameter("dim_Dm", unit(), self,
                "diffusivity of negative ions [m^2/s]")
        self.dim_Damb = daeParameter("dim_Damb", unit(), self,
                "ambipolar diffusivity [m^2/s]")
        self.Dp = daeParameter("Dp", unit(), self,
                "non-dimensional diffusivity of positive ions")
        self.Dm = daeParameter("Dm", unit(), self,
                "non-dimensional diffusivity of negative ions")
        self.rho_s = daeParameter("rho_s", unit(), self,
                "density of sites in solid [1/m^3]")
        self.csmax = daeParameter("csmax", unit(), self,
                "maximum concentration in solid [mol/m^3]")
        self.part_thickness = daeParameter("part_thickness",
                unit(), self,
                "Thickness of C3 particles")
        self.alpha = daeParameter("alpha", unit(), self,
                " Charge transfer coefficient")
        self.Tabs = daeParameter("Tabs", unit(), self,
                "Temperature in K")
        self.Tref = daeParameter("Tref", unit(), self,
                "Reference temperature in K")
        self.T = daeParameter("T", unit(), self,
                "Non dimensional temperature")
        self.k = daeParameter("k", unit(), self,
                "Boltzmann constant [J/(particle K)]")
        self.e = daeParameter("e", unit(), self,
                "Charge of proton, C")
        self.N_A = daeParameter("N_A", unit(), self,
                "Avogadro's number")
        self.F = daeParameter("F", unit(), self,
                "Faraday's number")
        self.C_rate = daeParameter("C_rate", unit(), self,
                "Discharge C-rate")
        self.currset = daeParameter("currset", unit(), self,
                "dimensionless current")
        self.dim_Vset = daeParameter("dim_Vset", unit(), self,
                "dimensional applied voltage (relative to " +
                "Delta V OCV of the  cell) [V]")
        self.Vset = daeParameter("Vset", unit(), self,
                "dimensionless applied voltage (relative to " +
                "Delta V OCV of the  cell)")
        self.cwet = daeParameter("c_wet", unit(), self,
                "Wetted surface concentration")
        self.dim_kappa = daeParameter("dim_kappa", unit(), self,
                "dimensional gradient penalty [J/m]")
        self.kappa = daeParameter("kappa", unit(), self,
                "kappa for each particle",
                [self.Ntrode, self.numpart])
        self.dim_a = daeParameter("dim_a", unit(), self,
                "dimensional reg sln [J]")
        self.a = daeParameter("a", unit(), self,
                "regular solution parameter for each particle [J]")
        self.dim_b = daeParameter("dim_b", unit(), self,
                "Stress coefficient [Pa]")
        self.b = daeParameter("b", unit(), self,
                "Stress coefficient for each particle")
        self.dim_k0 = daeParameter("dim_k0", unit(), self,
                "dimensional exchange current density rate constant [A/m^2]")
        self.k0 = daeParameter("k0", unit(), self,
                "exchange current density rate constant for each particle",
                [self.Ntrode, self.numpart])
        self.dim_mcond = daeParameter("dim_mcond", unit(), self,
                "dimensional conductivity of cathode [S/m]")
        self.mcond = daeParameter("mcond", unit(), self,
                "conductivity of cathode")
        self.dim_scond = daeParameter("dim_scond", unit(), self,
                "dimensional surface conductivity of particles [S]")
        self.scond = daeParameter("scond", unit(), self,
                "surface conductivity of particles",
                [self.Ntrode, self.numpart])
        self.Ltrode = daeParameter("Ltrode", unit(), self,
                "Length of the electrode")
        self.Lsep = daeParameter("Lsep", unit(), self,
                "Length of the separator")
        self.delx_sld = daeParameter("delx_sld", unit(), self,
                "size of discretization")
        self.Vstd = daeParameter("Vstd", unit(), self,
                "Standard potential [V]")
        self.psd_mean = daeParameter("psd_mean", unit(), self,
                "Particle size distribution mean [m]")
        self.psd_stddev = daeParameter("psd_stddev", unit(), self,
                "Particle size distribution stddev [m]")
        self.psd_num = daeParameter("psd_numVols", unit(), self,
                "Particle numbers of discretizations",
                [self.Ntrode, self.numpart])
        self.psd_len = daeParameter("psd_lengths", unit(), self,
                "Particle lengths [nm]",
                [self.Ntrode, self.numpart])
        self.psd_area = daeParameter("psd_active_areas", unit(), self,
                "Particle active areas [nm^2]",
                [self.Ntrode, self.numpart])
        self.psd_vol = daeParameter("psd_volumes", unit(), self,
                "Particle volumes [nm^3]",
                [self.Ntrode, self.numpart])
        self.type_ACR = daeParameter("type_ACR", unit(), self,
                "Bool: 1 for ACR type simulation")
        self.type_homog = daeParameter("type_homog", unit(), self,
                "Bool: 1 for homog type simulation")
        self.shape_sphere = daeParameter("shape_sphere", unit(), self,
                "Bool: 1 for spherical particles")
        self.shape_C3 = daeParameter("shape_C3", unit(), self,
                "Bool: 1 for C3 particles")
        self.cath_bulk_cond = daeParameter("cath_bulk_cond", unit(), self,
                "Bool: 1 to simulate bulk potential drop in cathode")
        self.cath_surf_cond = daeParameter("cath_surf_cond", unit(), self,
                "Bool: 1 to simulate surface potential drop in cathode")

    def DeclareEquations(self):
        daeModel.DeclareEquations(self)

        # Some values of domain lengths
        Ntrode = self.Ntrode.NumberOfPoints
        if Ntrode > 1:
            Nsep = self.Nsep.NumberOfPoints
        else:
            Nsep = 0
        Nlyte = Nsep + Ntrode
        numpart = self.numpart.NumberOfPoints
        Nsld_mat = np.zeros((Ntrode, numpart), dtype=np.integer)
        for i in range(Ntrode):
            for j in range(numpart):
                Nsld_mat[i, j] = self.Nsld_mat[i, j].NumberOfPoints

        # The porosity vector
        porosvec = np.empty(Nlyte + 1, dtype=object)
        # Use the Bruggeman relationship to approximate an effective
        # effect on the transport.
        porosvec[0:Nsep] = [self.poros_sep()**(3./2) for i in range(Nsep)]
        porosvec[Nsep:Nlyte+1] = [self.poros_trode()**(3./2) for i in
                range(Ntrode+1)]

#        # Prepare the noise
#        # maybe "numnoise" should be a parameter?
#        numnoise = tsteps/10
#        noise_prefac = 1e-3
#        noise_data = noise_prefac*np.random.randn(numnoise, Nsld)
#        # a vector going from 0 to the max simulation time.
#        time_vec = np.linspace(0, (1./self.currset.GetValue()), numnoise)
#        # Previous_output is common for all external functions
#        previous_output = []
#        # daeScalarExternalFunction (noise interpolation done as vector)
#        self.noise_local = np.empty(Nsld, dtype=object)
#        self.noise_local[:] = [noise("Noise", self, unit(), Time(),
#                                     time_vec, noise_data, previous_output, _position_)
#                               for _position_ in range(Nsld)]

        # Define the average concentration in each particle (algebraic
        # equations)
        for i in range(Ntrode):
            for j in range(numpart):
                eq = self.CreateEquation("cbar_vol{i}_part{j}".format(i=i,j=j))
                eq.Residual = (self.cbar_sld(i, j) -
                        Sum(self.c_sld[i, j].array([])) / Nsld_mat[i, j]
                        )
#                eq.BuildJacobianExpressions = True
                eq.CheckUnitsConsistency = False

        # Define the overall filling fraction in the cathode
        eq = self.CreateEquation("ffrac_cathode")
        eq.Residual = self.ffrac_cathode()
        numpartvol_tot = float(np.sum(Nsld_mat))
        for i in range(Ntrode):
            for j in range(numpart):
                eq.Residual -= (self.cbar_sld(i, j) *
                        (Nsld_mat[i, j]/numpartvol_tot))
        eq.CheckUnitsConsistency = False

        # Define dimensionless j_plus for each volume
        for i in range(Ntrode):
            eq = self.CreateEquation("j_plus_vol{i}".format(i=i))
            # Start with no reaction, then add reactions for each
            # particle in the volume.
            res = 0
            # sum over particle volumes in given electrode volume
            Vu = Sum(self.psd_vol.array(i, []))
            for  j in range(numpart):
                # The volume of this particular particle
                Vj = self.psd_vol(i, j)
                res += (Vj/Vu)*(Sum(self.c_sld[i, j].dt_array([])) /
                        Nsld_mat[i, j])
            eq.Residual = self.j_plus(i) - res
            eq.CheckUnitsConsistency = False

        # Calculate the solid concentration rates of change
        # (differential equations)
        for i in range(Ntrode):
            for j in range(numpart):
                # Prepare the RHS function
                Nij = Nsld_mat[i, j]
                RHS_c_sld_ij = self.calc_sld_dcs_dt(i, j)
                # Set up equations: dcdt = RHS
                for k in range(Nij):
                    eq = self.CreateEquation(
                            "dcsdt_vol{i}_part{j}_discr{k}".format(
                                i=i,j=j,k=k))
                    eq.Residual = self.c_sld[i, j].dt(k) - RHS_c_sld_ij[k]
                    eq.CheckUnitsConsistency = False

                # Also calculate the potential drop along cathode
                # particle surfaces, if desired
                simSurfCathCond = self.P.getboolean('Sim Params',
                        'simSurfCathCond')
                if simSurfCathCond:
                    # Conservation of charge in the solid particles with
                    # Ohm's Law
                    LHS = self.calc_part_surf_LHS(i, j)
                    k0_part = self.k0(i, j)
                    for k in range(Nij):
                        eq = self.CreateEquation(
                                "charge_cons_vol{i}_part{j}_discr{k}".format(
                                    i=i,j=j,k=k))
                        RHS = self.c_sld[i, j].dt(k) / k0_part
                        eq.Residual = LHS[k] - RHS
                        eq.CheckUnitsConsistency = False

        # Simulate the potential drop along the macroscopic-scale
        # cathode solid phase
        simBulkCathCond = self.P.getboolean('Sim Params',
                'simBulkCathCond')
        if simBulkCathCond:
            # Calculate the RHS for cathode conductivity
            phi_c = np.empty(Ntrode+2, dtype=object)
            phi_c[1:-1] = [self.phi_c(i) for i in range(Ntrode)]
            # No current passes into the electrolyte
            phi_c[0] = phi_c[1]
            # Potential at the current collector is set as a parameter
            phi_c[-1] = self.phi_cathode()
            dx = 1./Ntrode
            RHS_phi_c = -np.diff(-self.mcond()*np.diff(phi_c)/dx)/dx
        # Actually set up the equations for phi_c
        for i in range(Ntrode):
            eq = self.CreateEquation("phi_c{i}".format(i=i))
            if simBulkCathCond:
                eq.Residual = (-self.epsbeta()*self.j_plus(i) -
                        RHS_phi_c[i])
            else:
                eq.Residual = self.phi_c(i) - self.phi_cathode()

        # If we only have a single volume, electrolyte equations are
        # simple
        if Ntrode == 1:
            eq = self.CreateEquation("c_lyte")
            eq.Residual = self.c_lyte_trode.dt(0) - 0
            eq.CheckUnitsConsistency = False
            eq = self.CreateEquation("phi_lyte")
            eq.Residual = self.phi_lyte_trode(0) - self.phi_applied()
            eq.CheckUnitsConsistency = False
        else:
            # Calculate RHS for electrolyte equations
            Nlyte = Nsep + Ntrode
            c_lyte = np.empty(Nlyte, dtype=object)
            c_lyte[0:Nsep] = [self.c_lyte_sep(i) for i in range(Nsep)]
            c_lyte[Nsep:Nlyte] = [self.c_lyte_trode(i) for i in
                    range(Ntrode)]
            phi_lyte = np.empty(Nlyte, dtype=object)
            phi_lyte[0:Nsep] = [self.phi_lyte_sep(i) for i in range(Nsep)]
            phi_lyte[Nsep:Nlyte] = [self.phi_lyte_trode(i) for i in
                    range(Ntrode)]
            (RHS_c, RHS_phi) = self.calc_lyte_RHS(c_lyte, phi_lyte, Nlyte,
                    porosvec)

            # Equations governing the electrolyte in the separator
            for i in range(Nsep):
                # Mass Conservation
                eq = self.CreateEquation(
                        "sep_lyte_mass_cons_vol{i}".format(i=i))
                eq.Residual = (self.poros_sep()*self.c_lyte_sep.dt(i) -
                        RHS_c[i])
                eq.CheckUnitsConsistency = False
                # Charge Conservation
                eq = self.CreateEquation(
                        "sep_lyte_charge_cons_vol{i}".format(i=i))
                eq.Residual = (RHS_phi[i])
                eq.CheckUnitsConsistency = False
            # Equations governing the electrolyte in the electrode.
            # Here, we are coupled to the total reaction rates in the
            # solids.
            for i in range(Ntrode):
                # Mass Conservation
                eq = self.CreateEquation(
                        "trode_lyte_mass_cons_vol{i}".format(i=i))
                eq.Residual = (self.poros_trode()*self.c_lyte_trode.dt(i) +
                        self.epsbeta()*(1-self.tp())*self.j_plus(i) -
                        RHS_c[Nsep + i])
                eq.CheckUnitsConsistency = False
                # Charge Conservation
                eq = self.CreateEquation(
                        "trode_lyte_charge_cons_vol{i}".format(i=i))
                eq.Residual = (self.epsbeta()*self.j_plus(i) -
                        RHS_phi[Nsep + i])
                eq.CheckUnitsConsistency = False

        # Define the total current
        eq = self.CreateEquation("Total_Current")
        eq.Residual = self.current()
        dx = 1./Ntrode
        for i in range(Ntrode):
            eq.Residual -= dx*self.j_plus(i)
        eq.CheckUnitsConsistency = False

        if self.profileType == "CC":
            # Total Current Constraint Equation
            eq = self.CreateEquation("Total_Current_Constraint")
            eq.Residual = self.current() - self.currset()
            eq.CheckUnitsConsistency = False
        elif self.profileType == "CV":
            # Keep applied potential constant
            eq = self.CreateEquation("applied_potential")
            eq.Residual = self.phi_applied() - self.Vset()
            eq.CheckUnitsConsistency = False

    def calc_sld_dcs_dt(self, vol_indx, part_indx):
        # Get some useful information
        simSurfCathCond = self.P.getboolean('Sim Params',
                'simSurfCathCond')
        solidType = self.P.get('Sim Params', 'solidType')
        if simSurfCathCond and solidType != "ACR":
            raise Exception("Cannot do surface conductivity " +
                    "without ACR particles.")
        # shorthand
        i = vol_indx
        j = part_indx
        # Get variables for this particle/electrode volume
        phi_lyte = self.phi_lyte_trode(i)
        phi_m = self.phi_c(i)
        c_lyte = self.c_lyte_trode(i)
        # Get the relevant parameters for this particle
        k0 = self.k0(i, j)
        kappa = self.kappa(i, j) # only used for ACR
        cbar = self.cbar_sld(i, j) # only used for ACR
        # Number of volumes in current particle
        Nij = self.Nsld_mat[i, j].NumberOfPoints
        # Concentration (profile?) in the solid
        c_sld = np.empty(Nij, dtype=object)
        c_sld[:] = [self.c_sld[i, j](k) for k in range(Nij)]
#        # Potential (profile?) in the solid
#        if simSurfCathCond:
#            phi_m = np.empty(Nij
#        else:
#            phi_m = self.phi_c(i)
        # Calculate chemical potential of reduced state
        if solidType == "ACR":
            # Make a blank array to allow for boundary conditions
            cstmp = np.empty(Nij+2, dtype=object)
            cstmp[1:-1] = c_sld
            cstmp[0] = self.cwet()
            cstmp[-1] = self.cwet()
            dxs = 1./Nij
            curv = np.diff(cstmp, 2)/(dxs**2)
            mu_R = ( self.mu_reg_sln(c_sld) - kappa*curv
                    + self.b()*(c_sld - cbar) )
            # If we're also simulating potential drop along the solid,
            # use that instead of self.phi_c(i)
            if simSurfCathCond:
                phi_m = np.empty(Nij, dtype=object)
                phi_m[:] = [self.phi_sld[i, j](k) for k in range(Nij)]
        elif solidType == "homog":
            mu_R = self.mu_reg_sln(c_sld)
        act_R = np.exp(mu_R)
        gamma_ts = (1./(1-c_sld))
        # Assume dilute electrolyte
        act_O = c_lyte
        mu_O = np.log(act_O)
        # k0 is based on the _active_ area per volume for the region
        # of the solid of interest.
        ecd = ( k0 * act_O**(1-self.alpha())
                * act_R**(self.alpha()) / gamma_ts )
        # eta = electrochem pot_R - electrochem pot_O
        # eta = (mu_R + phi_R) - (mu_O + phi_O)
        eta = (mu_R + phi_m) - (mu_O + phi_lyte)
        return ( ecd*(np.exp(-self.alpha()*eta)
                - np.exp((1-self.alpha())*eta)) )

    def calc_lyte_RHS(self, cvec, phivec, Nlyte, porosvec):
        # The lengths are nondimensionalized by the electrode length
        dx = 1./self.Ntrode.NumberOfPoints
        # Mass conservation equations
        ctmp = np.empty(Nlyte + 2, dtype=object)
        ctmp[1:-1] = cvec
        # The total current flowing into the electrolyte is set
        ctmp[0] = (ctmp[1] +
                self.current()*self.epsbeta()*(1-self.tp())*dx
                )
        # No electrolyte flux at the separator
        ctmp[-1] = ctmp[-2]
        # Flux into the separator
        cflux = -porosvec*np.diff(ctmp)/dx
        # Divergence of the flux
        RHS_c = -np.diff(cflux)/dx

        # Charge conservation equations
        phitmp = np.empty(Nlyte + 2, dtype=object)
        phitmp[1:-1] = phivec
        # Currently, assume no rxn resistance at a lithium anode, and
        # measure relative to Li
        phitmp[0] = self.phi_applied()
        # No flux into cathode current collector from the electrolyte
        phitmp[-1] = phitmp[-2]
        # We need average values of c_lyte for the current densities
        # at the finite volume boundaries
        c_edges = (ctmp[0:-1] + ctmp[1:])/2.
        zp = self.zp()
        zm = self.zm()
        Dp = self.Dp()
        Dm = self.Dm()
        # Typo in Todd's code in currdens equation
        currdens = (-((Dp - Dm)*np.diff(ctmp)/dx) -
                (zp*Dp + zm*Dm)*c_edges*np.diff(phitmp)/dx)
        RHS_phi = -np.diff(porosvec*currdens)/dx
        return (RHS_c, RHS_phi)

    def calc_part_surf_LHS(self, vol_indx, part_indx):
        # shorthand
        i = vol_indx
        j = part_indx
        # Number of volumes in current particle
        Nij = self.Nsld_mat[i, j].NumberOfPoints
        # solid potential variables for this particle
        phi_tmp = np.empty(Nij + 2, dtype=object)
        phi_tmp[1:-1] = [self.phi_sld[i, j](k) for k in
                range(Nij)]
        # BC's -- "touching carbon at each end"
        phi_s_local = self.phi_c(i)
        phi_tmp[0] = phi_s_local
        phi_tmp[-1] = phi_s_local
        # LHS
        dx = 1./Nij
        phi_edges = (phi_tmp[0:-1] + phi_tmp[1:])/2.
#        curr_dens = -self.scond(i, j)*np.diff(phi_tmp, 1)/dx
        scond_vec = self.scond(i, j)*np.exp(-1*(phi_edges -
                phi_s_local))
        curr_dens = -scond_vec*np.diff(phi_tmp, 1)/dx
        return np.diff(curr_dens, 1)/dx

    def mu_reg_sln(self, c):
        return np.array([ self.a()*(1-2*c[i])
                + self.T()*Log(c[i]/(1-c[i]))
                for i in range(len(c)) ])

class simMPET(daeSimulation):
    def __init__(self, P=None):
        daeSimulation.__init__(self)
        if P is None:
            raise Exception("Need parameters input")
        self.P = P
        profileType = P.get('Sim Params', 'profileType')
        mean = P.getfloat('Sim Params', 'mean')
        stddev = P.getfloat('Sim Params', 'stddev')
        Ntrode = P.getint('Sim Params', 'Ntrode')
        numpart = P.getint('Sim Params', 'numpart')
        solidType = P.get('Sim Params', 'solidType')
        simSurfCathCond = P.getboolean('Sim Params', 'simSurfCathCond')
        # Make a length-sampled particle size distribution
        psd_raw = np.abs(stddev*np.random.randn(Ntrode, numpart) + mean)
        # For ACR particles, convert psd to integers -- number of steps
        if solidType == "ACR":
            solid_disc = P.getfloat('ACR info', 'solid_disc')
            self.psd_num = np.ceil(psd_raw/solid_disc).astype(np.integer)
            self.psd_len = solid_disc*self.psd_num
        # For homogeneous particles (only one "volume" per particle)
        elif solidType == "homog":
            # Each particle is only one volume
            self.psd_num = np.ones(psd_raw.shape).astype(np.integer)
            # The lengths are given by the original length distr.
            self.psd_len = psd_raw
        else:
            raise NotImplementedError("Input solidType not defined")
        # General parameters
        self.psd_mean = mean
        self.psd_stddev = stddev
        self.m = modMPET("mpet", Ntrode=Ntrode, numpart=numpart, P=P,
                simSurfCathCond=simSurfCathCond,
                profileType=profileType)

    def SetUpParametersAndDomains(self):
        # Extract info from the config file
        # Simulation
        dim_crate = self.P.getfloat('Sim Params', 'dim_crate')
        dim_Vset = self.P.getfloat('Sim Params', 'dim_Vset')
        Ntrode = self.P.getint('Sim Params', 'Ntrode')
        numpart = self.P.getint('Sim Params', 'numpart')
        T = self.P.getfloat('Sim Params', 'T')
        solidType = self.P.get('Sim Params', 'solidType')
        solidShape = self.P.get('Sim Params', 'solidShape')
        simBulkCathCond = self.P.getboolean('Sim Params',
                'simBulkCathCond')
        simSurfCathCond = self.P.getboolean('Sim Params',
                'simSurfCathCond')
        # Geometry
        Ltrode = self.P.getfloat('Geometry', 'Ltrode')
        Lsep = self.P.getfloat('Geometry', 'Lsep')
        Asep = self.P.getfloat('Geometry', 'Asep')
        Lp = self.P.getfloat('Geometry', 'Lp')
        poros = self.P.getfloat('Geometry', 'poros')
        # Electrolyte
        c0 = self.P.getfloat('Electrolyte Params', 'c0')
        zp = self.P.getfloat('Electrolyte Params', 'zp')
        zm = self.P.getfloat('Electrolyte Params', 'zm')
        Dp = self.P.getfloat('Electrolyte Params', 'Dp')
        Dm = self.P.getfloat('Electrolyte Params', 'Dm')
        # Cathode Material Properties
        dim_k0 = self.P.getfloat('Cathode Material Props', 'dim_k0')
        dim_a = self.P.getfloat('Cathode Material Props', 'dim_a')
        dim_kappa = self.P.getfloat('Cathode Material Props', 'dim_kappa')
        dim_b = self.P.getfloat('Cathode Material Props', 'dim_b')
        rhos = self.P.getfloat('Cathode Material Props', 'rhos')
        part_thick = self.P.getfloat('Cathode Material Props', 'part_thick')
        Vstd = self.P.getfloat('Cathode Material Props', 'Vstd')
        alpha = self.P.getfloat('Cathode Material Props', 'alpha')
        dim_mcond = self.P.getfloat('Cathode Material Props', 'dim_mcond')
        dim_scond = self.P.getfloat('Cathode Material Props', 'dim_scond')
        # ACR info
        cwet = self.P.getfloat('ACR info', 'cwet')
        solid_disc = self.P.getfloat('ACR info', 'solid_disc')
        # Constants
        k = self.P.getfloat('Constants', 'k')
        Tref = self.P.getfloat('Constants', 'Tref')
        e = self.P.getfloat('Constants', 'e')
        N_A = self.P.getfloat('Constants', 'N_A')
        # Calculated values
        # Faraday's number
        F = e*N_A
        # maximum concentration in cathode solid, mol/m^3
        csmax = rhos/N_A
        # Ambipolar diffusivity
        Damb = ((zp+zm)*Dp*Dm)/(zp*Dp+zm*Dm)
        # Cation transference number
        tp = zp*Dp / (zp*Dp + zm*Dm)
        # Diffusive time scale
        td = Ltrode**2 / Damb

        # Domains
        self.m.Ntrode.CreateArray(Ntrode)
        sep_frac = float(Lsep)/Ltrode
        Nsep = int(np.ceil(sep_frac*Ntrode))
        if Ntrode == 1:
            Nsep = 0
            sep_frac = 0
        else:
            sep_frac = float(Lsep)/Ltrode
            Nsep = int(np.ceil(sep_frac*Ntrode))
            self.m.Nsep.CreateArray(Nsep)
        self.m.numpart.CreateArray(numpart)
        for i in range(self.psd_num.shape[0]):
            for j in range(self.psd_num.shape[1]):
                self.m.Nsld_mat[i, j].CreateArray(int(self.psd_num[i, j]))

        # Parameters
        self.m.Tabs.SetValue(T)
        self.m.Tref.SetValue(Tref)
        self.m.T.SetValue(float(T)/Tref)
        self.m.k.SetValue(k)
        self.m.e.SetValue(e)
        self.m.N_A.SetValue(N_A)
        self.m.F.SetValue(F)
        self.m.alpha.SetValue(alpha)
        self.m.Ltrode.SetValue(Ltrode)
        self.m.Lsep.SetValue(Lsep)
        self.m.NumTrode.SetValue(Ntrode)
        self.m.NumSep.SetValue(Nsep)
        self.m.NumPart.SetValue(numpart)
        self.m.td.SetValue(td)
        self.m.zp.SetValue(zp)
        self.m.zm.SetValue(zm)
        self.m.tp.SetValue(tp)
        self.m.dim_Dp.SetValue(Dp)
        self.m.dim_Dm.SetValue(Dm)
        self.m.dim_Damb.SetValue(Damb)
        self.m.Dp.SetValue(Dp / Damb)
        self.m.Dm.SetValue(Dm / Damb)
        self.m.rho_s.SetValue(rhos)
        self.m.c_lyte0.SetValue(c0)
        self.m.csmax.SetValue(csmax)
        self.m.part_thickness.SetValue(part_thick)
        self.m.Lp.SetValue(Lp)
        self.m.dim_mcond.SetValue(dim_mcond)
        self.m.mcond.SetValue(dim_mcond * (td * k * N_A * Tref) /
                (Ltrode**2 * F**2 *c0))
        self.m.dim_scond.SetValue(dim_scond)
        self.m.poros_sep.SetValue(1.)
        self.m.poros_trode.SetValue(poros)
        self.m.epsbeta.SetValue((1-poros)*Lp*csmax/c0)
        self.m.phi_cathode.SetValue(0.)
        self.m.C_rate.SetValue(dim_crate)
        self.m.currset.SetValue(dim_crate*td/3600)
        self.m.dim_Vset.SetValue(dim_Vset)
        self.m.Vset.SetValue(dim_Vset*e/(k*Tref))
        self.m.dim_kappa.SetValue(dim_kappa)
        self.m.dim_k0.SetValue(dim_k0)
        self.m.dim_a.SetValue(dim_a)
        self.m.dim_b.SetValue(dim_b)
        self.m.a.SetValue(dim_a/(k*Tref))
        self.m.b.SetValue(dim_b/(k*Tref*rhos))
        self.m.psd_mean.SetValue(self.psd_mean)
        self.m.psd_stddev.SetValue(self.psd_stddev)
        if solidType == "ACR":
            self.m.type_ACR.SetValue(1.)
            self.m.type_homog.SetValue(0.)
        elif solidType == "homog":
            self.m.type_ACR.SetValue(0.)
            self.m.type_homog.SetValue(1.)
        if solidShape == "sphere":
            self.m.shape_sphere.SetValue(1.)
            self.m.shape_C3.SetValue(0.)
        if solidShape == "C3":
            self.m.shape_sphere.SetValue(0.)
            self.m.shape_C3.SetValue(1.)
        if simBulkCathCond:
            self.m.cath_bulk_cond.SetValue(1.)
        else:
            self.m.cath_bulk_cond.SetValue(0.)
        if simSurfCathCond:
            self.m.cath_surf_cond.SetValue(1.)
        else:
            self.m.cath_surf_cond.SetValue(0.)
        for i in range(Ntrode):
            for j in range(numpart):
                p_num = float(self.psd_num[i, j])
                p_len = self.psd_len[i, j]
                if solidShape == "sphere":
                    # Spherical particles
                    p_area = (4*np.pi)*p_len**2
                    p_vol = (4./3)*np.pi*p_len**3
                elif solidShape == "C3":
                    # C3 particles
                    p_area = 2 * 1.2263 * p_len**2
                    p_vol = 1.2263 * p_len**2 * part_thick
                else:
                    raise NotImplementedError("Input solidShape not defined")
                self.m.psd_num.SetValue(i, j, p_num)
                self.m.psd_len.SetValue(i, j, p_len)
                self.m.psd_area.SetValue(i, j, p_area)
                self.m.psd_vol.SetValue(i, j, p_vol)
                self.m.kappa.SetValue(i, j,
                        dim_kappa/(k*Tref*rhos*p_len**2))
                self.m.k0.SetValue(i, j,
                        ((p_area/p_vol)*dim_k0*td)/(F*csmax))
                self.m.scond.SetValue(i, j,
                        dim_scond * (k*Tref)/(dim_k0*e*p_len**2))
        self.m.cwet.SetValue(cwet)
        self.m.delx_sld.SetValue(solid_disc)
        self.m.Vstd.SetValue(Vstd)

    def SetUpVariables(self):
        Ntrode = self.m.Ntrode.NumberOfPoints
        if Ntrode > 1:
            Nsep = self.m.Nsep.NumberOfPoints
        else:
            Nsep = 0
        Nlyte = Nsep + Ntrode
        numpart = self.m.numpart.NumberOfPoints
        phi_cathode = self.m.phi_cathode.GetValue()
        # Set/guess values
        cs0 = self.P.getfloat('Sim Params', 'cs0')
        for i in range(Ntrode):
            # Guess initial volumetric reaction rates
            self.m.j_plus.SetInitialGuess(i, 0.0)
            # Guess initial value for the potential of the
            # cathode
            self.m.phi_c.SetInitialGuess(i, phi_cathode)
            for j in range(numpart):
                # Guess initial value for the average solid concentrations
                self.m.cbar_sld.SetInitialGuess(i, j, cs0)
                # Set initial solid concentration values
                Nij = self.m.Nsld_mat[i, j].NumberOfPoints
                for k in range(Nij):
                    self.m.c_sld[i, j].SetInitialCondition(k, cs0)
        # Set initial electrolyte concentration conditions
        c_lyte_init = 1.
        phi_guess = 0.
        for i in range(Nsep):
            self.m.c_lyte_sep.SetInitialCondition(i, c_lyte_init)
            self.m.phi_lyte_sep.SetInitialGuess(i, phi_guess)
        for i in range(Ntrode):
            self.m.c_lyte_trode.SetInitialCondition(i, c_lyte_init)
            self.m.phi_lyte_trode.SetInitialGuess(i, phi_guess)
        # Guess initial filling fraction
        self.m.ffrac_cathode.SetInitialGuess(cs0)
        # Guess the initial cell voltage
        self.m.phi_applied.SetInitialGuess(0.0)

class noise(daeScalarExternalFunction):
    def __init__(self, Name, Model, units, time, time_vec,
            noise_data, previous_output, position):
        arguments = {}
        self.counter = 0
        self.saved = 0
        self.previous_output = previous_output
        self.time_vec = time_vec
        self.noise_data = noise_data
        self.tlo = time_vec[0]
        self.thi = time_vec[-1]
        self.numnoise = len(time_vec)
        arguments["time"] = time
        self.position = position
        daeScalarExternalFunction.__init__(self, Name, Model, units, arguments)

    def Calculate(self, values):
        time = values["time"]
        # A derivative for Jacobian is requested - return always 0.0
        if time.Derivative != 0:
            return adouble(0)
        # Store the previous time value to prevent excessive
        # interpolation.
        if len(self.previous_output) > 0 and self.previous_output[0] == time.Value:
            self.saved += 1
            return adouble(float(self.previous_output[1][self.position]))
        indx = (float(time.Value - self.tlo)/(self.thi-self.tlo) *
                (self.numnoise - 1))
        ilo = np.floor(indx)
        ihi = np.ceil(indx)
        # If we're exactly at a time in time_vec
        if ilo == ihi:
            noise_vec = self.noise_data[ilo, :]
        else:
            noise_vec = (self.noise_data[ilo, :] +
                    (time.Value - self.time_vec[ilo]) /
                    (self.time_vec[ihi] - self.time_vec[ilo]) *
                    (self.noise_data[ihi, :] - self.noise_data[ilo, :])
                    )
        # previous_output is a reference to a common object and must
        # be updated here - not deleted.  using self.previous_output = []
        # it will delete the common object and create a new one
        self.previous_output[:] = [time.Value, noise_vec] # it is a list now not a tuple
        self.counter += 1
        return adouble(float(noise_vec[self.position]))

class MyMATDataReporter(daeMatlabMATFileDataReporter):
    """
    See Source code for pyDataReporting.daeMatlabMATFileDataReporter
    """
    def WriteDataToFile(self):
        mdict = {}
        for var in self.Process.Variables:
            mdict[var.Name] = var.Values
            mdict[var.Name + '_times'] = var.TimeValues
        try:
            scipy.io.savemat(self.ConnectionString,
                             mdict,
                             appendmat=False,
                             format='5',
                             long_field_names=False,
                             do_compression=False,
                             oned_as='row')
        except Exception, e:
            print 'Cannot call scipy.io.savemat(); is SciPy installed?\n' + str(e)

def setupDataReporters(simulation):
    """
    Create daeDelegateDataReporter and add data reporters:
     - daeMatlabMATFileDataReporter
    """
    datareporter = daeDelegateDataReporter()
    simulation.dr = MyMATDataReporter()
    datareporter.AddDataReporter(simulation.dr)
    # Connect data reporters
    simName = simulation.m.Name + strftime(" [%d.%m.%Y %H:%M:%S]", localtime())
    matfilename = os.path.join(os.getcwd(), "mpet_out.mat")
    if (simulation.dr.Connect(matfilename, simName) == False):
        sys.exit()
    return datareporter

def consoleRun(paramfile):
    # Prepare to interpret the config file for model inputs
    P = ConfigParser.RawConfigParser()
    P.read(paramfile)
    # Create Log, Solver, DataReporter and Simulation object
    log          = daePythonStdOutLog()
    daesolver    = daeIDAS()
    simulation   = simMPET(P)
    datareporter = setupDataReporters(simulation)

    # Use SuperLU direct sparse LA solver
    lasolver = pySuperLU.daeCreateSuperLUSolver()
#    lasolver = pyTrilinos.daeCreateTrilinosSolver("Amesos_Umfpack", "")
    daesolver.SetLASolver(lasolver)
    
    # Enable reporting of all variables
    simulation.m.SetReportingOn(True)

    # TODO -- optionally set the time horizon directly for CV?
    # Set the time horizon and the reporting interval
    # We need to get info about the system to figure out the
    # simulation time horizon
    dim_crate = P.getfloat('Sim Params', 'dim_crate')
    cs0 = P.getfloat('Sim Params', 'cs0')
    ffend = P.getfloat('Sim Params', 'ffend')
    tsteps = P.getfloat('Sim Params', 'tsteps')
    Ltrode = P.getfloat('Geometry', 'Ltrode')
    Dp = P.getfloat('Electrolyte Params', 'Dp')
    Dm = P.getfloat('Electrolyte Params', 'Dm')
    zp = P.getfloat('Electrolyte Params', 'zp')
    zm = P.getfloat('Electrolyte Params', 'zm')
    Damb = ((zp+zm)*Dp*Dm)/(zp*Dp+zm*Dm)
    td = Ltrode**2 / Damb
    currset = dim_crate * td/3600
    simulation.TimeHorizon = abs((ffend-cs0)/currset)
    simulation.ReportingInterval = simulation.TimeHorizon/tsteps

    # Connect data reporter
    simName = simulation.m.Name + strftime(" [%d.%m.%Y %H:%M:%S]", localtime())
    if(datareporter.Connect("", simName) == False):
        sys.exit()

    # Initialize the simulation
    simulation.Initialize(daesolver, datareporter, log)

    # Solve at time=0 (initialization)
    simulation.SolveInitial()

    # Run
    try:
        simulation.Run()
    except Exception as e:
        print str(e)
        simulation.ReportData(simulation.CurrentTime)
    simulation.Finalize()
    
if __name__ == "__main__":
    default_flag = 0
    default_file = "params_default.cfg"
    if len(sys.argv) < 2:
        default_flag = 1
        paramfile = default_file
    else:
        paramfile = sys.argv[1]
#    import timeit
#    time_tot = timeit.timeit("consoleRun()",
#            setup="from __main__ import consoleRun",  number=1)
#    print "Total time:", time_tot, "s"
#    import cProfile
#    cProfile.run("consoleRun()")
    consoleRun(paramfile)
    if default_flag:
        print "\n\n*** WARNING: Used default file, ""{fname}"" ***".format(
                fname=default_file)
        print "Pass other parameter file as an argument to this script\n"
    else:
        print "\n\nUsed parameter file ""{fname}""\n\n".format(
                fname=paramfile)
