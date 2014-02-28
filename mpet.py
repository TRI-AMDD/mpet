#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import sys
import os
import errno
import ConfigParser
import time

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

import mpet_params_IO
import delta_phi_fits

eps = 1e-12

# Define some variable types
mole_frac_t = daeVariableType(name="mole_frac_t", units=unit(),
        lowerBound=0, upperBound=1, initialGuess=0.25,
        absTolerance=1e-6)
elec_pot_t = daeVariableType(name="elec_pot_t", units=unit(),
        lowerBound=-1e20, upperBound=1e20, initialGuess=0,
        absTolerance=1e-5)

outdir_name = "sim_output"
outdir = os.path.join(os.getcwd(), outdir_name)

class modMPET(daeModel):
    def __init__(self, Name, Parent=None, Description="", D=None):
        daeModel.__init__(self, Name, Parent, Description)

        if (D is None):
            raise Exception("Need particle size distr. as input")
        self.D = D
        self.profileType = D['profileType']
        Nvol_c = D['Nvol_c']
        Nvol_s = D['Nvol_s']
        Nvol_a = D['Nvol_a']
        Npart_c = D['Npart_c']
        Npart_a = D['Npart_a']

        # Domains where variables are distributed
        if Nvol_s >= 1: # If we have a separator
            self.Nvol_s = daeDomain("Nvol_s", self, unit(),
                    "Number of control volumes in the separator")
#        if Nvol_a >= 1: # If we have an full anode
#            self.Nsep = daeDomain("Nsep", self, unit(),
#                    "Number of control volumes in the separator")
        self.Nvol_c = daeDomain("Nvol_c", self, unit(),
                "Number of control volumes in the electrode")
        self.Npart_c = daeDomain("Npart_c", self, unit(),
                "Number of particles sampled per electrode control volume")
        self.Nsld_mat = np.empty((Nvol_c, Npart_c), dtype=object)
        for i in range(Nvol_c):
            for j in range(Npart_c):
                self.Nsld_mat[i, j] = daeDomain("vol{i}_part{j}".format(
                    i=i, j=j), self, unit(),
                    "Number of discretizations for particle "
                    + "j in volume i".format(i=i,j=j))

        # Variables
        self.c_lyte_trode = daeVariable("c_lyte_trode", mole_frac_t, self,
                "Concentration in the electrolyte in the electrode",
                [self.Nvol_c])
        self.phi_lyte_trode = daeVariable("phi_lyte_trode", elec_pot_t, self,
                "Electrostatic potential in electrolyte in electrode",
                [self.Nvol_c])
        if Nvol_c > 1: # If we have a separator
            self.c_lyte_sep = daeVariable("c_lyte_sep", mole_frac_t, self,
                    "Concentration in the electrolyte in the separator",
                    [self.Nvol_s])
            self.phi_lyte_sep = daeVariable("phi_lyte_sep", elec_pot_t, self,
                    "Electrostatic potential in electrolyte in separator",
                    [self.Nvol_s])
        self.phi_applied = daeVariable("phi_applied", elec_pot_t, self,
                "Overall battery voltage (at anode current collector)")
        self.current = daeVariable("current", no_t, self,
                "Total current of the cell")
        self.c_sld = np.empty((Nvol_c, Npart_c), dtype=object)
        for i in range(Nvol_c):
            for j in range(Npart_c):
                self.c_sld[i, j] = daeVariable("solid_c_vol{i}_part{j}".format(
                    i=i, j=j), mole_frac_t, self,
                    "Concentration in each solid particle",
                    [self.Nsld_mat[i, j]])
        # Only make a variable of this if we have to -- it's a lot of
        # equations to keep track of for nothing if we don't need it.
        if D['simSurfCond_c']:
            self.phi_sld = np.empty((Nvol_c, Npart_c), dtype=object)
            for i in range(Nvol_c):
                for j in range(Npart_c):
                    self.phi_sld[i, j] = daeVariable("solid_p_vol{i}_part{j}".format(
                        i=i, j=j), elec_pot_t, self,
                        "Electrostatic potential in each solid particle",
                        [self.Nsld_mat[i, j]])
        self.cbar_sld = daeVariable("cbar_sld", mole_frac_t, self,
                "Average concentration in each particle",
                [self.Nvol_c, self.Npart_c])
        self.phi_c = daeVariable("phi_cath", elec_pot_t, self,
                "Electrostatic potential in the solid",
                [self.Nvol_c])
        self.j_plus = daeVariable("j_plus", no_t, self,
                "Rate of reaction of positives per solid volume",
                [self.Nvol_c])
        self.ffrac_cathode = daeVariable("ffrac_cathode",
                mole_frac_t, self,
                "Overall filling fraction of solids in cathode")

        # Parameters
        self.NumVol_c = daeParameter("NumVol_c", unit(), self,
                "Number of volumes in the electrode")
        self.NumPart_c = daeParameter("NumPart_c", unit(), self,
                "Number of particles in each electrode volume")
        self.NumVol_s = daeParameter("NumVol_s", unit(), self,
                "Number of volumes in the electrolyte")
        self.epsbeta = daeParameter("epsbeta", unit(), self,
                "porosity times beta")
        self.zp = daeParameter("zp", unit(), self,
                "cation charge number")
        self.zm = daeParameter("zm", unit(), self,
                "anion charge number")
        self.tp = daeParameter("tp", unit(), self,
                "positive transference number")
        self.poros_s = daeParameter("poros_s", unit(), self,
                "porosity in separator")
        self.poros_c = daeParameter("poros_c", unit(), self,
                "porosity in electrode")
        self.phi_cathode = daeParameter("phi_cathode", unit(), self,
                "potential at the cathode (phi_applied is relative to this)")
        self.td = daeParameter("td", unit(), self,
                "Diffusive time [s]")
        self.dim_Damb = daeParameter("dim_Damb", unit(), self,
                "ambipolar diffusivity [m^2/s]")
        self.dim_csmax_c = daeParameter("dim_csmax_c", unit(), self,
                "maximum lithium concentration in solid [mol/m^3]")
        self.Dp = daeParameter("Dp", unit(), self,
                "non-dimensional diffusivity of positive ions")
        self.Dm = daeParameter("Dm", unit(), self,
                "non-dimensional diffusivity of negative ions")
        self.Dsld_c = daeParameter("Dsld_c", unit(), self,
                "Diffusivity in cathode solid particles",
                [self.Nvol_c, self.Npart_c])
        self.alpha = daeParameter("alpha", unit(), self,
                " Charge transfer coefficient")
        self.T = daeParameter("T", unit(), self,
                "Non dimensional temperature")
        self.currset = daeParameter("currset", unit(), self,
                "dimensionless current")
        self.Vset = daeParameter("Vset", unit(), self,
                "dimensionless applied voltage (relative to " +
                "Delta V OCV of the  cell)")
        if self.D['etaFit_c']:
            self.dphi_eq_ref = daeParameter("dphi_eq_ref", unit(), self,
                    "dimensionless potential offset in referencing fit " +
                    "delta_phi_eq curves")
        self.cwet = daeParameter("c_wet", unit(), self,
                "Wetted surface concentration")
        self.kappa = daeParameter("kappa", unit(), self,
                "kappa for each particle",
                [self.Nvol_c, self.Npart_c])
        self.a = daeParameter("a", unit(), self,
                "regular solution parameter for each particle [J]",
                [self.Nvol_c, self.Npart_c])
        self.B = daeParameter("b", unit(), self,
                "Stress coefficient for each particle")
        self.k0 = daeParameter("k0", unit(), self,
                "exchange current density rate constant for each particle",
                [self.Nvol_c, self.Npart_c])
        self.lmbda = daeParameter("lmbda", unit(), self,
                "Marcus reorganizational energy")
        self.mcond = daeParameter("mcond", unit(), self,
                "conductivity of cathode")
        self.scond = daeParameter("scond", unit(), self,
                "surface conductivity of particles",
                [self.Nvol_c, self.Npart_c])
        self.psd_num = daeParameter("psd_numVols", unit(), self,
                "Particle numbers of discretizations",
                [self.Nvol_c, self.Npart_c])
        self.psd_len = daeParameter("psd_lengths", unit(), self,
                "Particle lengths [nm]",
                [self.Nvol_c, self.Npart_c])
        self.psd_area = daeParameter("psd_active_areas", unit(), self,
                "Particle active areas [nm^2]",
                [self.Nvol_c, self.Npart_c])
        self.psd_vol = daeParameter("psd_volumes", unit(), self,
                "Particle volumes [nm^3]",
                [self.Nvol_c, self.Npart_c])

    def DeclareEquations(self):
        daeModel.DeclareEquations(self)

        # Some values of domain lengths
        Nvol_c = self.Nvol_c.NumberOfPoints
        if Nvol_c > 1:
            Nvol_s = self.Nvol_s.NumberOfPoints
        else:
            Nvol_s = 0
        Nlyte = Nvol_s + Nvol_c
        Npart_c = self.Npart_c.NumberOfPoints
        Nsld_mat = np.zeros((Nvol_c, Npart_c), dtype=np.integer)
        for i in range(Nvol_c):
            for j in range(Npart_c):
                Nsld_mat[i, j] = self.Nsld_mat[i, j].NumberOfPoints

        # The porosity vector
        porosvec = np.empty(Nlyte + 1, dtype=object)
        # Use the Bruggeman relationship to approximate an effective
        # effect on the transport.
        porosvec[0:Nvol_s] = [self.poros_s()**(3./2) for i in range(Nvol_s)]
        porosvec[Nvol_s:Nlyte+1] = [self.poros_c()**(3./2) for i in
                range(Nvol_c+1)]

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
        for i in range(Nvol_c):
            for j in range(Npart_c):
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
        for i in range(Nvol_c):
            for j in range(Npart_c):
                eq.Residual -= (self.cbar_sld(i, j) *
                        (Nsld_mat[i, j]/numpartvol_tot))
        eq.CheckUnitsConsistency = False

        # Define dimensionless j_plus for each volume
        for i in range(Nvol_c):
            eq = self.CreateEquation("j_plus_vol{i}".format(i=i))
            # Start with no reaction, then add reactions for each
            # particle in the volume.
            res = 0
            # sum over particle volumes in given electrode volume
            Vu = Sum(self.psd_vol.array(i, []))
            for  j in range(Npart_c):
                # The volume of this particular particle
                Vj = self.psd_vol(i, j)
                res += (Vj/Vu)*(Sum(self.c_sld[i, j].dt_array([])) /
                        Nsld_mat[i, j])
            eq.Residual = self.j_plus(i) - res
            eq.CheckUnitsConsistency = False

        # Calculate the solid concentration rates of change
        # (differential equations)
        for i in range(Nvol_c):
            for j in range(Npart_c):
                # Prepare the RHS function
                Nij = Nsld_mat[i, j]
                (Mmat, RHS_c_sld_ij) = self.calc_sld_dcs_dt(i, j)
                dcdt_vec = np.empty(Nij, dtype=object)
                dcdt_vec[0:Nij] = [self.c_sld[i, j].dt(k) for k in range(Nij)]
                LHS_vec = self.MX(Mmat, dcdt_vec)
                # Set up equations: dcdt = RHS
                for k in range(Nij):
                    eq = self.CreateEquation(
                            "dcsdt_vol{i}_part{j}_discr{k}".format(
                                i=i,j=j,k=k))
#                    eq.Residual = self.c_sld[i, j].dt(k) - RHS_c_sld_ij[k]
                    eq.Residual = LHS_vec[k] - RHS_c_sld_ij[k]
                    eq.CheckUnitsConsistency = False

                # Also calculate the potential drop along cathode
                # particle surfaces, if desired
                simSurfCathCond = self.D['simSurfCond_c']
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
        simBulkCathCond = self.D['simBulkCond_c']
        if simBulkCathCond:
            # Calculate the RHS for cathode conductivity
            phi_c = np.empty(Nvol_c+2, dtype=object)
            phi_c[1:-1] = [self.phi_c(i) for i in range(Nvol_c)]
            # No current passes into the electrolyte
            phi_c[0] = phi_c[1]
            # Potential at the current collector is set as a parameter
            phi_c[-1] = self.phi_cathode()
            dx = 1./Nvol_c
            RHS_phi_c = -np.diff(-self.mcond()*np.diff(phi_c)/dx)/dx
        # Actually set up the equations for phi_c
        for i in range(Nvol_c):
            eq = self.CreateEquation("phi_c{i}".format(i=i))
            if simBulkCathCond:
                eq.Residual = (-self.epsbeta()*self.j_plus(i) -
                        RHS_phi_c[i])
            else:
                eq.Residual = self.phi_c(i) - self.phi_cathode()

        # If we only have a single volume, electrolyte equations are
        # simple
        if Nvol_c == 1:
            eq = self.CreateEquation("c_lyte")
            eq.Residual = self.c_lyte_trode.dt(0) - 0
            eq.CheckUnitsConsistency = False
            eq = self.CreateEquation("phi_lyte")
            eq.Residual = self.phi_lyte_trode(0) - self.phi_applied()
            eq.CheckUnitsConsistency = False
        else:
            # Calculate RHS for electrolyte equations
            Nlyte = Nvol_s + Nvol_c
            c_lyte = np.empty(Nlyte, dtype=object)
            c_lyte[0:Nvol_s] = [self.c_lyte_sep(i) for i in range(Nvol_s)]
            c_lyte[Nvol_s:Nlyte] = [self.c_lyte_trode(i) for i in
                    range(Nvol_c)]
            phi_lyte = np.empty(Nlyte, dtype=object)
            phi_lyte[0:Nvol_s] = [self.phi_lyte_sep(i) for i in range(Nvol_s)]
            phi_lyte[Nvol_s:Nlyte] = [self.phi_lyte_trode(i) for i in
                    range(Nvol_c)]
            (RHS_c, RHS_phi) = self.calc_lyte_RHS(c_lyte, phi_lyte, Nlyte,
                    porosvec)

            # Equations governing the electrolyte in the separator
            for i in range(Nvol_s):
                # Mass Conservation
                eq = self.CreateEquation(
                        "sep_lyte_mass_cons_vol{i}".format(i=i))
                eq.Residual = (self.poros_s()*self.c_lyte_sep.dt(i) -
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
            for i in range(Nvol_c):
                # Mass Conservation
                eq = self.CreateEquation(
                        "trode_lyte_mass_cons_vol{i}".format(i=i))
                eq.Residual = (self.poros_c()*self.c_lyte_trode.dt(i) +
                        self.epsbeta()*(1-self.tp())*self.j_plus(i) -
                        RHS_c[Nvol_s + i])
                eq.CheckUnitsConsistency = False
                # Charge Conservation
                eq = self.CreateEquation(
                        "trode_lyte_charge_cons_vol{i}".format(i=i))
                eq.Residual = (self.epsbeta()*self.j_plus(i) -
                        RHS_phi[Nvol_s + i])
                eq.CheckUnitsConsistency = False

        # Define the total current
        eq = self.CreateEquation("Total_Current")
        eq.Residual = self.current()
        dx = 1./Nvol_c
        for i in range(Nvol_c):
            eq.Residual -= dx*self.j_plus(i)
        eq.CheckUnitsConsistency = False

        if self.profileType == "CC":
            # Total Current Constraint Equation
            eq = self.CreateEquation("Total_Current_Constraint")
#            eq.Residual = self.current() - self.currset()
            timeHorizon = 1./Abs(self.currset())
            eq.Residual = self.current() - self.currset()*(1 -
                    np.exp(-Time()/(timeHorizon*1e-3)))
            eq.CheckUnitsConsistency = False
        elif self.profileType == "CV":
            # Keep applied potential constant
            eq = self.CreateEquation("applied_potential")
            eq.Residual = self.phi_applied() - self.Vset()
            eq.CheckUnitsConsistency = False

#        self.action = doNothingAction()
##        self.ON_CONDITION(Time() >= Constant(300*s),
#        self.ON_CONDITION(
##                Time() >= Constant(100*s) & Abs(self.phi_applied()) >= 60,
#                Abs(self.phi_applied()) >= 20,
#                switchToStates = [],
#                setVariableValues = [],
#                triggerEvents = [],
#                userDefinedActions = [self.action] )

    def calc_sld_dcs_dt(self, vol_indx, part_indx):
        # Get some useful information
        simSurfCathCond = self.D['simSurfCond_c']
        solidType_c = self.D['solidType_c']
        solidShape = self.D['solidShape_c']
        rxnType = self.D['rxnType_c']
        etaFit = self.D['etaFit_c']
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
        lmbda = self.lmbda() # Only used for Marcus
        alpha = self.alpha() # Only used for BV
        a = self.a(i, j)
        Ds = self.Dsld_c(i, j) # Only used for "diffn"
        # We need the (non-dimensional) temperature to get the
        # reaction rate dependence correct
        T = self.T()
        # Number of volumes in current particle
        Nij = self.Nsld_mat[i, j].NumberOfPoints
        # Concentration (profile?) in the solid
        c_sld = np.empty(Nij, dtype=object)
        c_sld[:] = [self.c_sld[i, j](k) for k in range(Nij)]
        # Calculate chemical potential of reduced state

        if solidType_c in ["ACR", "homog", "homog_sdn"]:
            if solidType_c == "ACR":
                # Make a blank array to allow for boundary conditions
                cstmp = np.empty(Nij+2, dtype=object)
                cstmp[1:-1] = c_sld
                cstmp[0] = self.cwet()
                cstmp[-1] = self.cwet()
                dxs = 1./Nij
                curv = np.diff(cstmp, 2)/(dxs**2)
                mu_R = ( self.mu_reg_sln(c_sld, a) - kappa*curv
                        + self.B()*(c_sld - cbar) )
                # If we're also simulating potential drop along the solid,
                # use that instead of self.phi_c(i)
                if simSurfCathCond:
                    phi_m = np.empty(Nij, dtype=object)
                    phi_m[:] = [self.phi_sld[i, j](k) for k in range(Nij)]
            elif solidType_c == "homog" or solidType_c == "homog_sdn":
                mu_R = self.mu_reg_sln(c_sld, a)
            # XXX -- Temp dependence!
            act_R = np.exp(mu_R)
            # Assume dilute electrolyte
            act_O = c_lyte
            mu_O = np.log(Max(eps, act_O))
            # eta = electrochem pot_R - electrochem pot_O
            # eta = (mu_R + phi_R) - (mu_O + phi_O)
            eta = (mu_R + phi_m) - (mu_O + phi_lyte)
            if rxnType == "Marcus":
                Rate = self.R_Marcus(k0, lmbda, c_lyte, c_sld, eta, T)
            elif rxnType == "BV":
                Rate = self.R_BV(k0, alpha, c_sld, act_O, act_R, eta, T)
            M = sprs.eye(Nij, format="csr")
            return (M, Rate)

        elif solidType_c in ["diffn"] and solidShape == "sphere":
            # For discretization background, see Zeng & Bazant 2013
            Rs = 1.
            dr = Rs/(Nij - 1)
            r_vec = np.linspace(0, Rs, Nij)
            vol_vec = r_vec**2 * dr + (1./12)*dr**3
            vol_vec[0] = (1./24)*dr**3
            vol_vec[-1] = (1./3)*(Rs**3 - (Rs - dr/2.)**3)
            M1 = sprs.diags([1./8, 3./4, 1./8], [-1, 0, 1],
                    shape=(Nij,Nij), format="csr")
            M1[1, 0] = M1[-2, -1] = 1./4
            M2 = sprs.diags(vol_vec, 0, format="csr")
            M = M1*M2
            RHS = np.empty(Nij, dtype=object)
            c_diffs = np.diff(c_sld)
            RHS[1:Nij - 1] = (
                    Ds*(r_vec[1:Nij - 1] + dr/2)**2*c_diffs[1:]/dr -
                    Ds*(r_vec[1:Nij - 1] - dr/2)**2*c_diffs[:-1]/dr )
            RHS[0] = Ds*(dr/2)**2*c_diffs[0]/dr
            # Figure out reaction rate information, assuming DILUTE
            # electrolyte AND solid
            # Take the surface concentration
            c_surf = c_sld[-1]
            # Overpotential
            delta_phi = phi_m - phi_lyte
            if etaFit:
                material = self.D['material_c']
                fits = delta_phi_fits.DPhiFits(self.D)
                phifunc = fits.materialData[material]
                delta_phi_eq = phifunc(c_surf, self.dphi_eq_ref())
            else:
                delta_phi_eq = T*np.log(Max(eps, c_lyte)/Max(eps, c_surf))
            eta = delta_phi - delta_phi_eq
            if rxnType == "Marcus":
                Rxn = self.R_Marcus(k0, lmbda, c_lyte, c_surf, eta, T)
            elif rxnType == "BV":
                # Assume dilute electrolyte
                act_O = c_lyte
                act_R = c_surf
                Rxn = self.R_BV(k0, alpha, c_surf, act_O, act_R, eta, T)
            RHS[-1] = (Rs**2 * Rxn -
                    Ds*(Rs - dr/2)**2*c_diffs[-1]/dr )
            return (M, RHS)
        elif solidType_c in ["diffn"] and solidShape == "C3":
            # TODO -- Implement
            raise
            return (M, RHS)

    def calc_lyte_RHS(self, cvec, phivec, Nlyte, porosvec):
        # The lengths are nondimensionalized by the electrode length
        dx = 1./self.Nvol_c.NumberOfPoints
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
        # XXX -- Temp dependence!
        scond_vec = self.scond(i, j)*np.exp(-1*(phi_edges -
                phi_s_local))
        curr_dens = -scond_vec*np.diff(phi_tmp, 1)/dx
        return np.diff(curr_dens, 1)/dx

    def mu_reg_sln(self, c, a):
        return np.array([ a*(1-2*c[i])
#                + self.T()*Log(c[i]/(1-c[i]))
#                + self.T()*Log((c[i]+eps)/(1-c[i]+eps))
                + self.T()*Log(Max(eps, c[i])/Max(eps, 1-c[i]))
                for i in range(len(c)) ])

    def R_BV(self, k0, alpha, c_sld, act_O, act_R, eta, T):
        gamma_ts = (1./(1-c_sld))
        ecd = ( k0 * act_O**(1-alpha)
                * act_R**(alpha) / gamma_ts )
        Rate = ( ecd *
            (np.exp(-alpha*eta/T) - np.exp((1-alpha)*eta/T)) )
        return Rate

    def R_Marcus(self, k0, lmbda, c_lyte, c_sld, eta, T):
        c_sld = np.array([Max(eps, c_sld[i]) for i in
            range(len(c_sld))])
#        alpha = 0.5*(1 + (T/lmbda) * np.log(Max(eps, c_lyte)/Max(eps, c_sld)))
        alpha = 0.5*(1 + (T/lmbda) * np.log(Max(eps, c_lyte)/c_sld))
        # We'll assume c_e = 1 (at the standard state for electrons)
#        ecd = ( k0 * np.exp(-lmbda/(4.*T)) *
        ecd = ( k0 *
                c_lyte**((3-2*alpha)/4.) *
                c_sld**((1+2*alpha)/4.) )
        Rate = ( ecd * np.exp(-eta**2/(4.*T*lmbda)) *
            (np.exp(-alpha*eta/T) - np.exp((1-alpha)*eta/T)) )
        return Rate

    def MX(self, mat, objvec):
        if type(mat) is not sprs.csr.csr_matrix:
            raise Exception("MX function designed for csr mult")
        n = objvec.shape[0]
        if (type(objvec[0]) == pyCore.adouble):
            out = np.empty(n, dtype=object)
        else:
            out = np.zeros(n, dtype=float)
        # Loop through the rows
        for i in range(n):
            low = mat.indptr[i]
            up = mat.indptr[i+1]
            if up > low:
                out[i] = np.sum(
                        mat.data[low:up] * objvec[mat.indices[low:up]] )
            else:
                out[i] = 0.0
        return out

class simMPET(daeSimulation):
    def __init__(self, D=None):
        daeSimulation.__init__(self)
        if D is None:
            raise Exception("Need parameters input")
        self.D = D
        self.test_input(D)
        mean_c = D['mean_c']
        stddev_c = D['stddev_c']
        Nvol_c = D['Nvol_c']
        Npart_c = D['Npart_c']
        solidType_c = D['solidType_c']
        # Make a length-sampled particle size distribution
        # Log-normally distributed
        if stddev_c == 0:
            psd_raw_c = mean*np.ones((Nvol_c, Npart_c))
        else:
            var = stddev_c**2
            mu = np.log((mean_c**2)/np.sqrt(var+mean_c**2))
            sigma = np.sqrt(np.log(var/(mean_c**2)+1))
            psd_raw_c = np.random.lognormal(mu, sigma,
                    size=(Nvol_c, Npart_c))
        # For particles with internal profiles, convert psd to
        # integers -- number of steps
        if solidType_c in ["ACR", "diffn"]:
            solidDisc_c = D['solidDisc_c']
            self.psd_num_c = np.ceil(psd_raw_c/solidDisc_c).astype(np.integer)
            self.psd_len_c = solidDisc_c*self.psd_num_c
        # For homogeneous particles (only one "volume" per particle)
        elif solidType_c in ["homog", "homog_sdn"]:
            # Each particle is only one volume
            self.psd_num_c = np.ones(psd_raw_c.shape).astype(np.integer)
            # The lengths are given by the original length distr.
            self.psd_len_c = psd_raw_c
        # General parameters
#        self.psd_mean = mean_c
#        self.psd_stddev = stddev_c
        self.m = modMPET("mpet", D=D)

    def SetUpParametersAndDomains(self):
        # Extract info from the config file
        # Simulation
        D = self.D
        Nvol_c = D['Nvol_c']
        Npart_c = D['Npart_c']
        solidType_c = D['solidType_c']
        solidShape_c = D['solidShape_c']
        # Geometry
        L_c = D['L_c']
        # Electrolyte
        zp = D['zp']
        zm = D['zm']
        Dp = D['Dp']
        Dm = D['Dm']
        # Cathode Material Properties
        # Cathode reaction
        # ACR info
        # Constants
        k = D['k']
        Tref = D['Tref']
        e = D['e']
        N_A = D['N_A']
        # Calculated values
        # Faraday's number
        F = e*N_A
        # maximum concentration in cathode solid, mol/m^3
        csmax_c = D['rhos_c']/N_A
        # Ambipolar diffusivity
        Damb = ((zp+zm)*Dp*Dm)/(zp*Dp+zm*Dm)
        # Cation transference number
        tp = zp*Dp / (zp*Dp + zm*Dm)
        # Diffusive time scale
        td = L_c**2 / Damb
        # Nondimensional Temperature
        T = float(D['Tabs'])/Tref

        # Domains
        self.m.Nvol_c.CreateArray(Nvol_c)
        sep_frac = float(D['L_s'])/L_c
        Nvol_s = int(np.ceil(sep_frac*Nvol_c))
        if Nvol_c == 1:
            Nvol_s = 0
            sep_frac = 0
        else:
            sep_frac = float(D['L_s'])/L_c
            Nvol_s = int(np.ceil(sep_frac*Nvol_c))
            self.m.Nvol_s.CreateArray(Nvol_s)
        self.m.Npart_c.CreateArray(Npart_c)
        for i in range(self.psd_num_c.shape[0]):
            for j in range(self.psd_num_c.shape[1]):
                self.m.Nsld_mat[i, j].CreateArray(int(self.psd_num_c[i, j]))

        # Parameters
        self.m.T.SetValue(T)
        self.m.alpha.SetValue(D['alpha_c'])
        self.m.NumVol_c.SetValue(Nvol_c)
        self.m.NumVol_s.SetValue(Nvol_s)
        self.m.NumPart_c.SetValue(Npart_c)
        self.m.td.SetValue(td)
        self.m.zp.SetValue(zp)
        self.m.zm.SetValue(zm)
        self.m.tp.SetValue(tp)
        self.m.dim_csmax_c.SetValue(csmax_c)
        self.m.dim_Damb.SetValue(Damb)
        self.m.Dp.SetValue(Dp / Damb)
        self.m.Dm.SetValue(Dm / Damb)
        self.m.mcond.SetValue(D['mcond_c'] * (td * k * N_A * Tref) /
                (L_c**2 * F**2 *D['c0']))
        self.m.poros_s.SetValue(1.)
        self.m.poros_c.SetValue(D['poros_c'])
        self.m.epsbeta.SetValue((1-D['poros_c'])*D['PL_c']*csmax_c/D['c0'])
        self.m.phi_cathode.SetValue(0.)
        self.m.currset.SetValue(D['Crate']*td/3600)
        self.m.Vset.SetValue(D['Vset']*e/(k*Tref))
        if self.D['etaFit_c']:
            material = self.D['material_c']
            fits = delta_phi_fits.DPhiFits(self.D)
            phifunc = fits.materialData[material]
            self.m.dphi_eq_ref.SetValue(phifunc(self.D['cs0_c'], 0))
        self.m.lmbda.SetValue(D['lambda_c']/(k*Tref))
        self.m.B.SetValue(D['B_c']/(k*Tref*D['rhos_c']))
        for i in range(Nvol_c):
            for j in range(Npart_c):
                p_num = float(self.psd_num_c[i, j])
                p_len = self.psd_len_c[i, j]
                # k0 is based on the _active_ area per volume for the region
                # of the solid of interest.
                if solidShape_c == "sphere":
                    # Spherical particles
                    p_area = (4*np.pi)*p_len**2
                    p_vol = (4./3)*np.pi*p_len**3
                elif solidShape_c == "C3":
                    # C3 particles
                    p_area = 2 * 1.2263 * p_len**2
                    p_vol = 1.2263 * p_len**2 * D['partThick_c']
                self.m.psd_num.SetValue(i, j, p_num)
                self.m.psd_len.SetValue(i, j, p_len)
                self.m.psd_area.SetValue(i, j, p_area)
                self.m.psd_vol.SetValue(i, j, p_vol)
                self.m.kappa.SetValue(i, j,
                        D['kappa_c']/(k*Tref*D['rhos_c']*p_len**2))
                self.m.k0.SetValue(i, j,
                        ((p_area/p_vol)*D['k0_c']*td)/(F*csmax_c))
                self.m.scond.SetValue(i, j,
                        D['scond_c'] * (k*Tref)/(D['k0_c']*e*p_len**2))
                self.m.Dsld_c.SetValue(i, j,
                        D['Dsld_c']*(p_area/p_vol)*td/p_len)
                if solidType_c in ["homog", "ACR",  "diffn"]:
                    self.m.a.SetValue(i, j, D['Omega_a_c']/(k*Tref))
                elif solidType_c == "homog_sdn":
                    # Not sure about factor of nondimensional T. Thus,
                    # only use this when T = 1, Tabs = Tref = 298
                    self.m.a.SetValue(i, j, T*self.size2regsln(p_len))
        self.m.cwet.SetValue(D['cwet_c'])

    def SetUpVariables(self):
        Nvol_c = self.m.Nvol_c.NumberOfPoints
        if Nvol_c > 1:
            Nvol_s = self.m.Nvol_s.NumberOfPoints
        else:
            Nvol_s = 0
        Nlyte = Nvol_s + Nvol_c
        Npart_c = self.m.Npart_c.NumberOfPoints
        phi_cathode = self.m.phi_cathode.GetValue()
        # Set/guess values
        cs0_c = self.D['cs0_c']
        for i in range(Nvol_c):
            # Guess initial volumetric reaction rates
            self.m.j_plus.SetInitialGuess(i, 0.0)
            # Guess initial value for the potential of the
            # cathode
            self.m.phi_c.SetInitialGuess(i, phi_cathode)
            for j in range(Npart_c):
                # Guess initial value for the average solid concentrations
                self.m.cbar_sld.SetInitialGuess(i, j, cs0_c)
                # Set initial solid concentration values
                Nij = self.m.Nsld_mat[i, j].NumberOfPoints
                for k in range(Nij):
                    self.m.c_sld[i, j].SetInitialCondition(k, cs0_c)
        # Set initial electrolyte concentration conditions
        c_lyte_init = 1.
        phi_guess = 0.
        for i in range(Nvol_s):
            self.m.c_lyte_sep.SetInitialCondition(i, c_lyte_init)
            self.m.phi_lyte_sep.SetInitialGuess(i, phi_guess)
        for i in range(Nvol_c):
            self.m.c_lyte_trode.SetInitialCondition(i, c_lyte_init)
            self.m.phi_lyte_trode.SetInitialGuess(i, phi_guess)
        # Guess initial filling fraction
        self.m.ffrac_cathode.SetInitialGuess(cs0_c)
        # Guess the initial cell voltage
        self.m.phi_applied.SetInitialGuess(0.0)

    def size2regsln(self, size):
        """
        This function returns the non-dimensional regular solution
        parameter which creates a barrier height that corresponds to
        the given particle size (C3 particle, measured in nm in the
        [100] direction). The barrier height vs size is taken from
        Cogswell 2013, and the reg sln vs barrier height was done by
        TRF 2014.
        """
        # First, this function wants the argument to be in [nm]
        size *= 1e+9
        # Parameters for polynomial curve fit
        p1 = -1.168e4
        p2 = 2985
        p3 = -208.3
        p4 = -8.491
        p5 = -10.25
        p6 = 4.516
        # The nucleation barrier depends on the ratio of the particle
        # wetted area to total particle volume.
        # *Wetted* area to volume ratio for C3 particles (Cogswell
        # 2013 or Kyle Smith)
        AV = 3.6338/size
        # Fit function (TRF, "SWCS" paper 2014)
        param = p1*AV**5 + p2*AV**4 + p3*AV**3 + p4*AV**2 + p5*AV + p6
        if param < 2:
            param = 2
#        param = [param[i] if param[i] >= 2 else 2 for i in
#                range(len(param))]
        return param

    def test_input(self, D):
        solidType_c = D['solidType_c']
        solidShape_c = D['solidShape_c']
        if D['Tabs'] != 298:
            raise Exception("Temp dependence not implemented")
        if D['simSurfCond_c'] and solidType_c != "ACR":
            raise Exception("simSurfCond req. ACR")
        if solidType_c in ["ACR", "homog_sdn"] and solidShape_c != "C3":
            raise Exception("ACR and homog_sdn req. C3 shape")
        if solidType_c not in ["ACR", "homog", "homog_sdn", "diffn"]:
            raise NotImplementedError("Input solidType_c not defined")
        if solidShape_c not in ["C3", "sphere"]:
            raise NotImplementedError("Input solidShape_c not defined")
        if solidType_c == "homog_sdn" and (D['Tabs'] != 298 or
                D['Tref'] != 298):
            raise NotImplementedError("homog_snd req. Tref=Tabs=298")
        if solidType_c in ["diffn"] and solidShape_c != "sphere":
            raise NotImplementedError("diffn currently req. sphere")
        if D['etaFit_c'] and solidType_c != "diffn":
            raise NotImplementedError("etafit req. solidType_c = diffn")
        return

#    def Run(self):
#        """
#        Overload the simulation "run" function so that the simulation
#        terminates when the specified condition is satisfied.
#        """
#        while self.CurrentTime < self.TimeHorizon:
#            t_step = self.CurrentTime + self.ReportingInterval
#            if t_step > self.TimeHorizon:
#                t_step = self.TimeHorizon
#
#            self.Log.Message("Integrating from %.2f to %.2fs ..." % (self.CurrentTime, t_step), 0)
#            self.IntegrateUntilTime(t_step, eStopAtModelDiscontinuity)
#            self.ReportData(self.CurrentTime)
#
#            if self.LastSatisfiedCondition:
#                self.Log.Message('Condition: [{0}] satisfied at time {1}s'.format(self.LastSatisfiedCondition, self.CurrentTime), 0)
#                self.Log.Message('Stopping the simulation...', 0)
#                return

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

class doNothingAction(daeAction):
    def __init__(self):
        daeAction.__init__(self)
    def Execute(self):
        pass

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
    simName = simulation.m.Name + time.strftime(" [%d.%m.%Y %H:%M:%S]",
            time.localtime())
    matDataName = "output_data.mat"
    matfilename = os.path.join(outdir, matDataName)
    if (simulation.dr.Connect(matfilename, simName) == False):
        sys.exit()
    return datareporter

def consoleRun(D):
    # Create Log, Solver, DataReporter and Simulation object
    log          = daePythonStdOutLog()
    daesolver    = daeIDAS()
    simulation   = simMPET(D)
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
    Dp = D['Dp']
    Dm = D['Dm']
    zp = D['zp']
    zm = D['zm']
    Damb = ((zp+zm)*Dp*Dm)/(zp*Dp+zm*Dm)
    td = D['L_c']**2 / Damb
    currset = D['Crate'] * td/3600.
    if D['profileType'] == "CC":
        simulation.TimeHorizon = D['capFrac']/currset
    else: # CV simulation
        simulation.TimeHorizon = D['tend']/td
    simulation.ReportingInterval = simulation.TimeHorizon/D['tsteps']

    # Connect data reporter
    simName = simulation.m.Name + time.strftime(" [%d.%m.%Y %H:%M:%S]",
            time.localtime())
    if(datareporter.Connect("", simName) == False):
        sys.exit()

    # Initialize the simulation
    simulation.Initialize(daesolver, datareporter, log)

#    # Save model report
#    simulation.m.SaveModelReport(simulation.m.Name + ".xml")

    # Solve at time=0 (initialization)
    simulation.SolveInitial()

    # Run
    try:
        simulation.Run()
#    except Exception as e:
    except Exception as e:
        print str(e)
        simulation.ReportData(simulation.CurrentTime)
    except KeyboardInterrupt:
        print "\nphi_applied at ctrl-C:", simulation.m.phi_applied.GetValue(), "\n"
        simulation.ReportData(simulation.CurrentTime)
    simulation.Finalize()

if __name__ == "__main__":
    timeStart = time.time()
    default_flag = 0
    default_file = "params_default.cfg"
    if len(sys.argv) < 2:
        default_flag = 1
        paramfile = default_file
    else:
        paramfile = sys.argv[1]
    # Get the parameters dictionary (and the config instance) from the
    # parameter file
    IO = mpet_params_IO.mpetIO()
    P = IO.getConfig(paramfile)
    D = IO.getDictFromConfig(P)
    # Make sure there's a place to store the output
    try:
        os.makedirs(outdir)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
    paramFileName = "output_params.cfg"
    paramFile = os.path.join(outdir, paramFileName)
    IO.writeConfigFile(P, filename=paramFile)
    consoleRun(D)
    if default_flag:
        print "\n\n*** WARNING: Used default file, ""{fname}"" ***".format(
                fname=default_file)
        print "Pass other parameter file as an argument to this script\n"
    else:
        print "\n\nUsed parameter file ""{fname}""\n\n".format(
                fname=paramfile)
    timeEnd = time.time()
    print "Total time:", (timeEnd - timeStart), "s"
