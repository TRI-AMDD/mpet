#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import sys
import os
from time import localtime, strftime
#import math

import numpy as np
import scipy.sparse as sprs
import scipy.interpolate as sint

from daetools.pyDAE import *
from daetools.pyDAE.data_reporters import *
from daetools.solvers.superlu import pySuperLU
#from daetools.solvers.superlu_mt import pySuperLU_MT
from daetools.solvers.trilinos import pyTrilinos
#from daetools.solvers.intel_pardiso import pyIntelPardiso
from pyUnits import s
#from pyUnits import s, kg, m, K, Pa, mol, J, W

#########################################################################
# CONSTANTS
k = 1.381e-23      # Boltzmann constant
T = 298            # Temp, K
Tref = 298         # Reference temp, K (for non-dimensionalization)
e = 1.602e-19      # Charge of proton, C
N_A = 6.02e23      # Avogadro's number
F = e*N_A          # Faraday's number

#########################################################################
# SET DIMENSIONAL VALUES HERE
#init_voltage = 3.5                 # Initial cell voltage, V
dim_crate = 5                   # Battery discharge c-rate
#currset = 0.001                         # Battery discharge c-rate

# Electrode properties
Ltrode = 50e-6                      # electrode thickness, m
Lsep = 25e-6                        # separator thickness, m
Asep = 1e-4                         # area of separator, m^2
Lp = 0.69                           # Volume loading percent active material
poros = 0.4                         # Porosity
c0 = 1000                           # Initial electrolyte conc., mol/m^3
zp = 1                              # Cation charge number
zm = 1                              # Anion charge number
Dp = 2.2e-10                        # Cation diff, m^2/s, LiPF6 in EC/DMC
Dm = 2.94e-10                       # Anion diff, m^2/s, LiPF6 in EC/DMC
Damb = ((zp+zm)*Dp*Dm)/(zp*Dp+zm*Dm)    # Ambipolar diffusivity
tp = zp*Dp / (zp*Dp + zm*Dm)        # Cation transference number

# Particle size distribution
mean = 20e-9                       # Average particle size, m
stddev = 2e-9                      # Standard deviation, m

# Material properties
dim_k0 = 0.1                       # Exchange current density, A/m^2 (0.1 for h2/Pt)
dim_a = 1.8560e-20                 # Regular solution parameter, J
dim_kappa = 5.0148e-10             # Gradient penalty, J/m
dim_b = 0.1916e9                   # Stress, Pa
rhos = 1.3793e28                   # site density, 1/m^3
csmax = rhos/N_A                   # maximum concentration, mol/m^3
cwet = 0.98                        # Dimensionless wetted conc.
#wet_thick = 2e-9                   # Thickness of wetting on surf.
part_thick = 20e-9                 # C3 Particle thickness, m
Vstd = 3.422                       # Standard potential, V
alpha = 0.5                        # Charge transfer coefficient

# Discretization settings
Ntrode = 6                              # Numer disc. in x direction (volumes in electrode)
solid_disc = 1e-9                   # Discretization size of solid, m (MUST BE LESS THAN LAMBDA)
numpart = 4                         # Particles per volume
tsteps = 200                        # Number disc. in time


# Define some variable types
mole_frac_t = daeVariableType(name="mole_frac_t", units=unit(),
        lowerBound=0, upperBound=1, initialGuess=0.25,
        absTolerance=1e-6)
elec_pot_t = daeVariableType(name="elec_pot_t", units=unit(),
        lowerBound=-1e20, upperBound=1e20, initialGuess=0,
        absTolerance=1e-5)

class modMPET(daeModel):
    def __init__(self, Name, Parent=None, Description="", psd=None):
        print "Init modMPET"
        daeModel.__init__(self, Name, Parent, Description)

        # TODO -- change psd to simply Num vols and numpart
        if psd is None:
            raise Exception("Need particle size distr. as input")

        # Domains where variables are distributed
        self.Nsep = daeDomain("Nsep", self, unit(),
                "Number of control volumes in the separator")
        self.Ntrode = daeDomain("Ntrode", self, unit(),
                "Number of control volumes in the electrode")
        self.numpart = daeDomain("numpart", self, unit(),
                "Number of particles sampled per electrode control volume")
#        self.Nsld_max = daeDomain("Nsld_max", self, unit(),
#                "Maximum number of discretizations in a particle")
#        self.Nsld = daeDomain("solids", self, unit(),
#                "Total number of stored solid concentration values")
        self.Nsld_mat = np.empty(psd.shape, dtype=object)
        for i in range(psd.shape[0]):
            for j in range(psd.shape[1]):
                self.Nsld_mat[i, j] = daeDomain("vol{i}_part{j}".format(
                    i=i, j=j), self, unit(),
                    "Number of discretizations for particle "
                    + "j in volume i".format(i=i,j=j))

        # Variables
        self.c_lyte_sep = daeVariable("c_lyte_sep", mole_frac_t, self,
                "Concentration in the electrolyte in the separator",
                [self.Nsep])
        self.c_lyte_trode = daeVariable("c_lyte_trode", mole_frac_t, self,
                "Concentration in the electrolyte in the electrode",
                [self.Ntrode])
        self.phi_lyte_sep = daeVariable("phi_lyte_sep", elec_pot_t, self,
                "Electrostatic potential in electrolyte in separator",
                [self.Nsep])
        self.phi_lyte_trode = daeVariable("phi_lyte_trode", elec_pot_t, self,
                "Electrostatic potential in electrolyte in electrode",
                [self.Ntrode])
        self.phi_applied = daeVariable("phi_applied", elec_pot_t, self,
                "Overall battery voltage (at anode current collector)")
        self.c_sld = np.empty(psd.shape, dtype=object)
        for i in range(psd.shape[0]):
            for j in range(psd.shape[1]):
                self.c_sld[i, j] = daeVariable("solid_vol{i}_part{j}".format(
                    i=i, j=j), mole_frac_t, self,
                    "Concentration in the solid particle",
                    [self.Nsld_mat[i, j]])
        self.cbar_sld = daeVariable("cbar_sld", mole_frac_t, self,
                "Average concentration in each particle",
                [self.Ntrode, self.numpart])
#        self.phi_sld = daeVariable("phi_sld", elec_pot_t, self,
#                "Electrostatic potential in the solid")
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
#                [self.Ntrode, self.numpart])
        self.dim_b = daeParameter("dim_b", unit(), self,
                "Stress coefficient [Pa]")
        self.b = daeParameter("b", unit(), self,
                "Stress coefficient for each particle")
#                [self.Ntrode, self.numpart])
#        self.k0 = daeParameter("k0", unit(), self,
#                "exchange current density rate constant")
        self.dim_k0 = daeParameter("dim_k0", unit(), self,
                "dimensional exchange current density rate constant [A/m^2]")
        self.k0 = daeParameter("k0", unit(), self,
                "exchange current density rate constant for each particle",
                [self.Ntrode, self.numpart])
#        self.aO = daeParameter("aO", unit(), self,
#                "activity of the oxidized state")
#        self.Lx = daeParameter("Lx", unit(), self,
#                "Length of particle")
#        self.Ntrode = daeParameter("Ntrode", unit(), self,
#                "Number of volumes in particle")
        self.Ltrode = daeParameter("Ltrode", unit(), self,
                "Length of the electrode")
        self.Lsep = daeParameter("Lsep", unit(), self,
                "Length of the separator")
        self.delx_sld = daeParameter("delx_sld", unit(), self,
                "size of discretization")
#        self.part_sizes = daeParameter("sizes", unit(), self,
#                "particle sizes, in discretization counts",
#                [self.Nx, self.numpart])
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

    def DeclareEquations(self):
        print "DeclareEquations()"
        daeModel.DeclareEquations(self)

        # Some values of domain lengths
        Nsep = self.Nsep.NumberOfPoints
        Ntrode = self.Ntrode.NumberOfPoints
        Nlyte = Nsep + Ntrode
        numpart = self.numpart.NumberOfPoints
#        print type(self.numpart)
#        print numpart
#        print type(numpart)
#        Nsld = self.Nsld.NumberOfPoints
#        Nsld_max = self.Nsld_max.NumberOfPoints
        Nsld_mat = np.zeros((Ntrode, numpart), dtype=np.integer)
        for i in range(Ntrode):
            for j in range(numpart):
#                print type(self.Nsld_mat[i, j])
                Nsld_mat[i, j] = self.Nsld_mat[i, j].NumberOfPoints
#                print Nsld_mat[i, j]
#                print type(Nsld_mat[i, j])

        # The porosity vector
#        porosvec = np.empty(Nlyte, dtype=object)
#        # Use the Bruggeman relationship to approximate an effective
#        # effect on the transport.
#        porosvec[0:Nsep] = [self.poros_sep()**(3./2) for i in range(Nsep)]
#        porosvec[Nsep:Nlyte] = [self.poros_trode()**(3./2) for i in
#                range(Ntrode)]
        porosvec = np.empty(Nlyte + 1, dtype=object)
        # Use the Bruggeman relationship to approximate an effective
        # effect on the transport.
        porosvec[0:Nsep] = [self.poros_sep()**(3./2) for i in range(Nsep)]
        porosvec[Nsep:Nlyte+1] = [self.poros_trode()**(3./2) for i in
                range(Ntrode+1)]

#        # Prepare the noise
#        # XXX -- maybe this should be a parameter?
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
#                print eq.Residual

        # Define the overall filling fraction in the cathode
        # XXX -- This is wrong: need to adjust for particle sizes
        eq = self.CreateEquation("ffrac_cathode")
        eq.Residual = (self.ffrac_cathode() -
                Sum(self.cbar_sld.array([], []))/(Ntrode*numpart)
                )
        eq.CheckUnitsConsistency = False
#        print eq.Residual

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
#                res += (Vj/Vu)*self.cbar_sld.dt(i, j)
                res += (Vj/Vu)*(Sum(self.c_sld[i, j].dt_array([])) /
                        Nsld_mat[i, j])
            eq.Residual = self.j_plus(i) - res
#            # DEBUG -- single particle per volume
#            j = 0
#            eq.Residual = (self.j_plus(i) -
#                    Sum(self.c_sld[i, j].dt_array([])) /
#                        Nsld_mat[i, j])
            eq.CheckUnitsConsistency = False
#            print eq.Residual

        # Calculate the solid concentration rates of change
        # (differential equations)
        for i in range(Ntrode):
            phi_lyte = self.phi_lyte_trode(i)
            c_lyte = self.c_lyte_trode(i)
            for j in range(numpart):
                # Prepare the RHS function
                Nij = Nsld_mat[i, j]
                c_sld = np.empty(Nij, dtype=object)
                c_sld[:] = [self.c_sld[i, j](k) for k in range(Nij)]
#                k0 = self.k0(i, j)
#                kappa = self.kappa(i, j)
                # TODO -- have dcd_dt functions extract all needed
                # information -- no arguments
                # For ACR style particles
                RHS_c_sld_ij = self.calc_ACR_dcs_dt(c_sld, phi_lyte,
                        c_lyte, i, j)
#                # For homogeneous particles
#                k0 = self.k0(i, j)
#                RHS_c_sld_ij = self.calc_homog_dcs_dt(c_sld, phi_lyte,
#                        c_lyte, k0)
                # Set up equations: dcdt = RHS
                for k in range(Nij):
                    eq = self.CreateEquation(
                            "dcsdt_vol{i}part{j}_{k}".format(
                                i=i,j=j,k=k))
                    eq.Residual = self.c_sld[i, j].dt(k) - RHS_c_sld_ij[k]
                    eq.CheckUnitsConsistency = False
#                    print eq.Residual

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
#            print eq.Residual
            # Charge Conservation
            eq = self.CreateEquation(
                    "sep_lyte_charge_cons_vol{i}".format(i=i))
            eq.Residual = (RHS_phi[i])
            eq.CheckUnitsConsistency = False
#            print eq.Residual
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
#            eq.Residual = self.poros_trode()*self.c_lyte_trode.dt(i)
#            Vu = Sum(self.psd_vol.array(i, []))
#            for j in range(numpart):
#                Vj = self.psd_vol(i, j)
#                eq.Residual += (self.epsbeta()*(1-self.tp()) * (Vj/Vu) *
#                        Sum(self.c_sld[i, j].dt_array([]))/Nsld_mat[i, j])
#            eq.Residual -= RHS_c[Nsep + i]
            eq.CheckUnitsConsistency = False
#            print eq.Residual
            # Charge Conservation
            eq = self.CreateEquation(
                    "trode_lyte_charge_cons_vol{i}".format(i=i))
            eq.Residual = (self.epsbeta()*self.j_plus(i) -
                    RHS_phi[Nsep + i])
#            eq.Residual = -RHS_phi[Nsep + i]
#            for j in range(numpart):
#                Vj = self.psd_vol(i, j)
#                eq.Residual += (self.epsbeta() * (Vj/Vu) *
#                        Sum(self.c_sld[i, j].dt_array([]))/Nsld_mat[i, j])
            eq.CheckUnitsConsistency = False
#            print eq.Residual

        # Total Current Constraint Equation
        eq = self.CreateEquation("Total_Current_Constraint")
#        eq.Residual = (
#                Sum(self.j_plus.array([])) - self.currset())
        eq.Residual = -self.currset()
        dx = 1./Ntrode
        for i in range(Ntrode):
            eq.Residual += dx*self.j_plus(i)
#        eq.Residual = -self.currset()
#        for i in range(Ntrode):
#            Vu = Sum(self.psd_vol.array(i, []))
#            for j in range(numpart):
#                Vj = self.psd_vol(i, j)
#                eq.Residual += ((Vj/Vu) *
#                        Sum(self.c_sld[i, j].dt_array([]))/Nsld_mat[i, j])
        eq.CheckUnitsConsistency = False
#        print eq.Residual

    def calc_ACR_dcs_dt(self, cs, phi_lyte, c_lyte, vol_indx, part_indx):
#        part_steps = self.N.NumberOfPoints
        # shorthand
        i = vol_indx
        j = part_indx
        # Get the relevant parameters for this particle
        k0 = self.k0(i, j)
        kappa = self.kappa(i, j)
        cbar = self.cbar_sld(i, j)
        # Number of volumes in current particle
        part_steps = len(cs)
        # Make a blank array to allow for boundary conditions
        cstmp = np.empty(part_steps+2, dtype=object)
        cstmp[1:-1] = cs
        cstmp[0] = self.cwet()
        cstmp[-1] = self.cwet()
        dxs = 1./part_steps
        curv = np.diff(cstmp, 2)/(dxs**2)
        mu = ( self.mu_reg_sln(cs) - kappa*curv
                + self.b()*(cs - cbar) )
        act_R = np.exp(mu)
        gamma_ts = (1./(1-cs))
        # Assume dilute electrolyte
        act_O = c_lyte
        # k0 is based on the _active_ area per volume for the region
        # of the solid of interest.
        ecd = ( k0 * act_O**(1-self.alpha())
                * act_R**(self.alpha()) / gamma_ts )
#        delta_phi = phi_sld - phi_lyte
        delta_phi = self.phi_cathode() - phi_lyte
        eta = mu + delta_phi
        return ( ecd*(np.exp(-self.alpha()*eta)
                - np.exp((1-self.alpha())*eta)) )

    def calc_homog_dcs_dt(self, c_sld, phi_lyte, c_lyte, k0):
        mu = self.mu_reg_sln(c_sld)
        act_R = np.exp(mu)
        gamma_ts = (1./(1-c_sld))
        act_O = c_lyte
        ecd = ( k0 * act_O**(1-self.alpha())
                * act_R**(self.alpha()) / gamma_ts )
        delta_phi = self.phi_cathode() - phi_lyte
        eta = mu + delta_phi
        return ( ecd*(np.exp(-self.alpha()*eta) -
                np.exp((1-self.alpha())*eta)) )

    def calc_lyte_RHS(self, cvec, phivec, Nlyte, porosvec):
#        Nlyte = Nsep + Ntrode
#        print porosvec
        # The lengths are nondimensionalized by the electrode length
        dx = 1./Ntrode
        # Mass conservation equations
#        ctmp = np.empty(Nlyte + 1, dtype=object)
#        ctmp[:-1] = cvec
        ctmp = np.empty(Nlyte + 2, dtype=object)
        ctmp[1:-1] = cvec
        # XXX -- Is this overspecifying the system? No.
        # The total current flowing into the electrolyte is set
        ctmp[0] = (ctmp[1] +
                self.currset()*self.epsbeta()*(1-self.tp())*dx
                )
#        ctmp[0] = ctmp[1]
        # No electrolyte flux at the separator
        ctmp[-1] = ctmp[-2]
        # Flux into the separator
        cflux = -porosvec*np.diff(ctmp)/dx
        # Divergence of the flux
        RHS_c = -np.diff(cflux)/dx
#        RHS_c = 0*RHS_c

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

    def mu_reg_sln(self, c):
        if (type(c[0]) == pyCore.adouble):
            isAdouble = True
        else:
            isAdouble = False
        if isAdouble:
            return np.array([ self.a()*(1-2*c[i])
                    + self.T()*Log(c[i]/(1-c[i]))
                    for i in range(len(c)) ])
        else:
            return ( self.a.GetValue()*(1-2*c)
                    + self.T.GetValue()*np.log(c/(1-c)) )

class simMPET(daeSimulation):
    def __init__(self, psd, psd_mean, psd_stddev):
        print "initalize simulation"
        daeSimulation.__init__(self)
        if psd is None:
            raise Exception("Need particle size distr. as input")
        self.psd = psd
        self.psd_mean = psd_mean
        self.psd_stddev = psd_stddev
        self.m = modMPET("mpet", psd=psd)

    def SetUpParametersAndDomains(self):
        print "SetUpParametersAndDomains"
        # Domains
        Ntrode = self.psd.shape[0]
        numpart = self.psd.shape[1]
        self.m.Ntrode.CreateArray(Ntrode)
        sep_frac = float(Lsep)/Ltrode
        Nsep = int(np.ceil(sep_frac*Ntrode))
        self.m.Nsep.CreateArray(Nsep)
        self.m.numpart.CreateArray(numpart)
        for i in range(self.psd.shape[0]):
            for j in range(self.psd.shape[1]):
                self.m.Nsld_mat[i, j].CreateArray(int(self.psd[i, j]))

        # Parameters
        self.m.Tabs.SetValue(T)
        self.m.Tref.SetValue(Tref)
        self.m.T.SetValue(float(T)/Tref)
        self.m.k.SetValue(k)
        self.m.e.SetValue(e)
        self.m.N_A.SetValue(N_A)
        self.m.F.SetValue(F)
        self.m.alpha.SetValue(alpha)
        td = Ltrode**2 / Damb
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
        csmax = rhos/N_A
        self.m.c_lyte0.SetValue(c0)
        self.m.csmax.SetValue(csmax)
        self.m.part_thickness.SetValue(part_thick)
        self.m.Lp.SetValue(Lp)
        self.m.poros_sep.SetValue(1.)
        self.m.poros_trode.SetValue(poros)
        self.m.epsbeta.SetValue((1-poros)*Lp*csmax/c0)
        self.m.phi_cathode.SetValue(0.)
        self.m.C_rate.SetValue(dim_crate)
        self.m.currset.SetValue(dim_crate*td/3600)
        self.m.dim_kappa.SetValue(dim_kappa)
        self.m.dim_k0.SetValue(dim_k0)
        self.m.dim_a.SetValue(dim_a)
        self.m.dim_b.SetValue(dim_b)
        self.m.a.SetValue(dim_a/(k*Tref))
        self.m.b.SetValue(dim_b/(k*Tref*rhos))
        self.m.psd_mean.SetValue(self.psd_mean)
        self.m.psd_stddev.SetValue(self.psd_stddev)
        for i in range(Ntrode):
            for j in range(numpart):
                p_num = float(self.psd[i, j])
                p_len = p_num*solid_disc
#                # Spherical particles
#                p_area = (4*np.pi)*p_len**2
#                p_vol = (4./3)*np.pi*p_len**3
                # C3 particles
                p_area = 2 * 1.2263 * p_len**2
                p_vol = 1.2263 * p_len**2 * part_thick
                self.m.psd_num.SetValue(i, j, p_num)
                self.m.psd_len.SetValue(i, j, p_len)
                self.m.psd_area.SetValue(i, j, p_area)
                self.m.psd_vol.SetValue(i, j, p_vol)
                self.m.kappa.SetValue(i, j,
                        dim_kappa/(k*Tref*rhos*p_len**2))
                self.m.k0.SetValue(i, j,
                        ((p_area/p_vol)*dim_k0*td)/(F*csmax))
        self.m.cwet.SetValue(cwet)
        self.m.delx_sld.SetValue(solid_disc)
        self.m.Vstd.SetValue(Vstd)

    def SetUpVariables(self):
        print "SetUpVariables"
        Nsep = self.m.Nsep.NumberOfPoints
        Ntrode = self.m.Ntrode.NumberOfPoints
        Nlyte = Nsep + Ntrode
        numpart = self.m.numpart.NumberOfPoints
        # Set/guess values
        cs0 = 0.01
        for i in range(Ntrode):
            # Guess initial volumetric reaction rates
            self.m.j_plus.SetInitialGuess(i, 0.0)
            for j in range(numpart):
                # Set initial value for the average solid concentrations
#                self.m.cbar_sld.SetInitialCondition(i, j, cs0)
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
        self.m.ffrac_cathode.SetInitialGuess(0.01)
        # Guess the initial cell voltage
        self.m.phi_applied.SetInitialGuess(0.0)
        print "done SetUpVariables"

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
#        print 'At time = %f the function with index = %d interpolated the noise' % (time.Value, int(self.position)) 
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
    print "setupDataReporter()"
    datareporter = daeDelegateDataReporter()
    simulation.dr = MyMATDataReporter()
    datareporter.AddDataReporter(simulation.dr)
    # Connect data reporters
    simName = simulation.m.Name + strftime(" [%d.%m.%Y %H:%M:%S]", localtime())
    matfilename = os.path.join(os.getcwd(),
            "acr_sp_{curr}C.mat".format(curr=dim_crate))
#    print matfilename
    if (simulation.dr.Connect(matfilename, simName) == False):
        sys.exit()
    return datareporter

def consoleRun():
    print "START"
    # A bit awkward, but we'll have to generate our
    # particle size distribution here and pass it to the
    # simulation/model.
    psd_raw = stddev*np.abs(np.random.randn(Ntrode, numpart)) + mean
    # Convert psd to integers -- number of steps
    psd = np.ceil(psd_raw/solid_disc).astype(np.integer)
#    # For homogeneous particles
#    psd = (50*np.ones((Ntrode, numpart))).astype(np.integer)
    # TODO -- pass mean and stddev into simulation/model and store as
    # parameters
    # Create Log, Solver, DataReporter and Simulation object
    log          = daePythonStdOutLog()
    daesolver    = daeIDAS()
    simulation   = simMPET(psd, mean, stddev)
    datareporter = setupDataReporters(simulation)

    # Use SuperLU direct sparse LA solver
    print "Set up LU solver"
    lasolver = pySuperLU.daeCreateSuperLUSolver()
#    lasolver = pyTrilinos.daeCreateTrilinosSolver("Amesos_Umfpack", "")
    daesolver.SetLASolver(lasolver)
    
    # Enable reporting of all variables
    simulation.m.SetReportingOn(True)

    # Set the time horizon and the reporting interval
    print "set up simulation time parameters"
#    simulation.TimeHorizon = 0.98/abs(dim_crate)
    td = Ltrode**2 / Damb
    currset = dim_crate*td/3600
#    simulation.TimeHorizon = 1./abs(currset)
    simulation.TimeHorizon = 0.95/abs(currset)
#    simulation.TimeHorizon = 0.98/abs(0.001)
#    simulation.TimeHorizon = 0.50/abs(0.001)
    simulation.ReportingInterval = simulation.TimeHorizon/tsteps

    # Connect data reporter
    print "connect to data reporter"
    simName = simulation.m.Name + strftime(" [%d.%m.%Y %H:%M:%S]", localtime())
    if(datareporter.Connect("", simName) == False):
        sys.exit()

    # Initialize the simulation
    print "initialize"
    simulation.Initialize(daesolver, datareporter, log)

    # Solve at time=0 (initialization)
    print "Initialize the system at t = 0"
    simulation.SolveInitial()

    # Run
    print "run the simulation"
    try:
        simulation.Run()
    except Exception as e:
        print str(e)
        simulation.ReportData(simulation.CurrentTime)
    simulation.Finalize()
    
#    print 'Number of numpy.interp1d calls in noise ext. functions:'
#    print simulation.m.noise_local.counter
#    print simulation.m.noise_local.saved
#    print [n.counter for n in simulation.m.noise_vec]

if __name__ == "__main__":
    import timeit
    time_tot = timeit.timeit("consoleRun()",
            setup="from __main__ import consoleRun",  number=1)
    print "Total time:", time_tot, "s"
#    import cProfile
#    cProfile.run("consoleRun()")
#    consoleRun()
