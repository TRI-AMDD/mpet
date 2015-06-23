#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import sys
import os
import shutil
import errno
import ConfigParser
import time
import subprocess

import numpy as np
import scipy.sparse as sprs
import scipy.special as spcl
import scipy.interpolate as sint
#import scipy.io as sio

from daetools.pyDAE import *
from daetools.pyDAE.data_reporters import *
from daetools.solvers.superlu import pySuperLU
#from daetools.solvers.superlu_mt import pySuperLU_MT
#from daetools.solvers.trilinos import pyTrilinos
#from daetools.solvers.intel_pardiso import pyIntelPardiso
#from pyUnits import s
#from pyUnits import s, kg, m, K, Pa, mol, J, W

import mpet_params_IO
import delta_phi_fits

eps = -1e-12

# Define some variable types
mole_frac_t = daeVariableType(name="mole_frac_t", units=unit(),
        lowerBound=0, upperBound=1, initialGuess=0.25,
        absTolerance=1e-7)
elec_pot_t = daeVariableType(name="elec_pot_t", units=unit(),
        lowerBound=-1e20, upperBound=1e20, initialGuess=0,
        absTolerance=1e-5)


class modMPET(daeModel):
    def __init__(self, Name, Parent=None, Description="", D=None):
        daeModel.__init__(self, Name, Parent, Description)

        if (D is None):
            raise Exception("Need particle size distr. as input")
        self.D = D
        self.profileType = D['profileType']
        Nvol_s = D['Nvol_s']
        Nvol_ac = np.array([D['Nvol_ac'][0], D['Nvol_ac'][1]])
        Npart_ac = np.array([D['Npart_ac'][0], D['Npart_ac'][1]])
        if Nvol_ac[0] >= 1: # If we have a full anode
            trodes = [0, 1]
        else:
            trodes = [1]
        self.trodes = trodes

        # Domains where variables are distributed
        self.Ntrode = daeDomain("Ntrode", self, unit(),
                "Number of porous electrodes in the cell")
        if Nvol_s >= 1: # If we have a separator
            self.Nvol_s = daeDomain("Nvol_s", self, unit(),
                    "Number of control volumes in the separator")
        self.Nvol_ac = np.empty(2, dtype=object)
        self.Npart_ac = np.empty(2, dtype=object)
        self.Nsld_mat_ac = np.empty(2, dtype=object)
        for l in trodes:
            self.Nvol_ac[l] = daeDomain("Nvol_{l}".format(l=l),
                    self, unit(),
                    "Number of control volumes in the electrode " +
                    "{l}".format(l=l))
            self.Npart_ac[l] = daeDomain("Npart_{l}".format(l=l),
                    self, unit(),
                    "Number of particles sampled per control " +
                    "volume in electrode {l}".format(l=l))
            Nvol = Nvol_ac[l]
            Npart = Npart_ac[l]
            Nsld_mat = np.empty((Nvol, Npart), dtype=object)
            for i in range(Nvol):
                for j in range(Npart):
                    Nsld_mat[i, j] = daeDomain("trode{l}_vol{i}_part{j}".format(
                        i=i, j=j, l=l), self, unit(),
                        "Number of discretizations for particle "
                        + "j in volume i".format(i=i,j=j))
            self.Nsld_mat_ac[l] = Nsld_mat

        # Variables
        self.c_lyte_ac = np.empty(2, dtype=object)
        self.phi_lyte_ac = np.empty(2, dtype=object)
        self.c1_sld_ac = np.empty(2, dtype=object)
        self.c2_sld_ac = np.empty(2, dtype=object)
        self.phi_sld_ac = np.empty(2, dtype=object)
        self.c1bar_sld_ac = np.empty(2, dtype=object)
        self.c2bar_sld_ac = np.empty(2, dtype=object)
        self.phi_ac = np.empty(2, dtype=object)
        self.phipart_ac = np.empty(2, dtype=object)
        self.j_plus_ac = np.empty(2, dtype=object)
        self.ffrac_ac = np.empty(2, dtype=object)
        for l in trodes:
            # Concentration/potential in electrode regions of elyte
            self.c_lyte_ac[l] = daeVariable("c_lyte_{l}".format(l=l),
                    mole_frac_t, self,
                    "Concentration in the electrolyte in " +
                    "electrode {l}".format(l=l),
                    [self.Nvol_ac[l]])
            self.phi_lyte_ac[l] = daeVariable("phi_lyte_{l}".format(l=l),
                    elec_pot_t, self,
                    "Electrostatic potential in electrolyte in " +
                    "electrode {l}".format(l=l),
                    [self.Nvol_ac[l]])
            # Concentration in electrode active particles
            Nvol = Nvol_ac[l]
            Npart = Npart_ac[l]
            self.c1_sld_ac[l] = np.empty((Nvol, Npart), dtype=object)
            self.c2_sld_ac[l] = np.empty((Nvol, Npart), dtype=object)
            for i in range(Nvol):
                for j in range(Npart):
                    self.c1_sld_ac[l][i, j] = daeVariable(
                            "c1_sld_trode{l}vol{i}part{j}".format(
                            i=i, j=j, l=l), mole_frac_t, self,
                            "Concentration in each solid particle",
                            [self.Nsld_mat_ac[l][i, j]])
                    self.c2_sld_ac[l][i, j] = daeVariable(
                            "c2_sld_trode{l}vol{i}part{j}".format(
                            i=i, j=j, l=l), mole_frac_t, self,
                            "Concentration in each solid particle",
                            [self.Nsld_mat_ac[l][i, j]])
            # Potential in electrode active particles
            # Only make a variable of solid potentials if we have to
            # -- it's a lot of equations to keep track of for nothing
            # if we don't need it.
            if D['simSurfCond_ac'][l]:
                self.phi_sld_ac[l] = np.empty((Nvol, Npart), dtype=object)
                for i in range(Nvol):
                    for j in range(Npart):
                        self.phi_sld_ac[l][i, j] = daeVariable(
                                "p_sld_trode{l}_vol{i}_part{j}".format(
                                i=i, j=j, l=l), elec_pot_t, self,
                                "Electrostatic potential in each solid particle",
                                [self.Nsld_mat[l][i, j]])
            else:
                self.phi_sld_ac[l] = False
            # Average active particle concentrations
            self.c1bar_sld_ac[l] = daeVariable("c1bar_sld_{l}".format(l=l),
                    mole_frac_t, self,
                    "Average concentration in each particle",
                    [self.Nvol_ac[l], self.Npart_ac[l]])
            self.c2bar_sld_ac[l] = daeVariable("c2bar_sld_{l}".format(l=l),
                    mole_frac_t, self,
                    "Average concentration in each particle",
                    [self.Nvol_ac[l], self.Npart_ac[l]])
            self.phi_ac[l] = daeVariable("phi_{l}".format(l=l),
                    elec_pot_t, self,
                    "Electrostatic potential in the bulk solid",
                    [self.Nvol_ac[l]])
            self.phipart_ac[l] = daeVariable("phipart_{l}".format(l=l),
                    elec_pot_t, self,
                    "Electrostatic potential at each particle",
                    [self.Nvol_ac[l], self.Npart_ac[l]])
            self.j_plus_ac[l] = daeVariable("j_plus_{l}".format(l=l),
                    no_t, self,
                    "Rate of reaction of positives per solid volume",
                    [self.Nvol_ac[l]])
            self.ffrac_ac[l] = daeVariable("ffrac_{l}".format(l=l),
                mole_frac_t, self,
                "Overall filling fraction of solids in electrodes")
        if Nvol_s >= 1: # If we have a separator
            self.c_lyte_s = daeVariable("c_lyte_s", mole_frac_t, self,
                    "Concentration in the electrolyte in the separator",
                    [self.Nvol_s])
            self.phi_lyte_s = daeVariable("phi_lyte_s", elec_pot_t, self,
                    "Electrostatic potential in electrolyte in separator",
                    [self.Nvol_s])
        self.phi_applied = daeVariable("phi_applied", elec_pot_t, self,
                "Overall battery voltage (at anode current collector)")
        self.current = daeVariable("current", no_t, self,
                "Total current of the cell")

        # Parameters
        self.NumTrode = daeParameter("NumTrode", unit(), self,
                "Number of electroes simulated (1 or 0)")
        self.NumVol_ac = daeParameter("NumVol_ac", unit(), self,
                "Number of volumes in the electrode",
                [self.Ntrode])
        self.NumPart_ac = daeParameter("NumPart_ac", unit(), self,
                "Number of particles in each electrode volume",
                [self.Ntrode])
        self.NumVol_s = daeParameter("NumVol_s", unit(), self,
                "Number of volumes in the electrolyte")
        self.L_ac = daeParameter("L_ac", unit(), self,
                "Length of electrodes (ndim to L_c)",
                [self.Ntrode])
        self.L_s = daeParameter("L_s", unit(), self,
                "Length of separator (ndim to L_c)")
        self.epsbeta_ac = daeParameter("epsbeta_ac", unit(), self,
                "porosity times beta in electrodes",
                [self.Ntrode])
        self.zp = daeParameter("zp", unit(), self,
                "cation charge number")
        self.zm = daeParameter("zm", unit(), self,
                "anion charge number")
        self.tp = daeParameter("tp", unit(), self,
                "positive transference number")
        self.poros_s = daeParameter("poros_s", unit(), self,
                "porosity in separator")
        self.poros_ac = daeParameter("poros_ac", unit(), self,
                "porosity in electrode",
                [self.Ntrode])
        self.phi_cathode = daeParameter("phi_cathode", unit(), self,
                "reference potential, at the cathode + "
                "(phi_applied is relative to this)")
        self.td = daeParameter("td", unit(), self,
                "Diffusive time [s]")
        self.dim_Damb = daeParameter("dim_Damb", unit(), self,
                "ambipolar diffusivity [m^2/s]")
        self.Dp = daeParameter("Dp", unit(), self,
                "non-dimensional diffusivity of positive ions")
        self.Dm = daeParameter("Dm", unit(), self,
                "non-dimensional diffusivity of negative ions")
        self.Dsld_ac = np.empty(2, dtype=object)
        self.kappa_ac = np.empty(2, dtype=object)
        self.Omga_ac = np.empty(2, dtype=object)
        self.Omgb_ac = np.empty(2, dtype=object)
        self.Omgc_ac = np.empty(2, dtype=object)
        self.k0_ac = np.empty(2, dtype=object)
        self.beta_s_ac = np.empty(2, dtype=object)
        self.delta_L_ac = np.empty(2, dtype=object)
        self.MHC_Aa_ac = np.empty(2, dtype=object)
        self.scond_ac = np.empty(2, dtype=object)
        self.G_ac = np.empty(2, dtype=object)
        self.psd_num_ac = np.empty(2, dtype=object)
        self.psd_len_ac = np.empty(2, dtype=object)
        self.psd_area_ac = np.empty(2, dtype=object)
        self.psd_vol_ac = np.empty(2, dtype=object)
        for l in trodes:
            self.Dsld_ac[l] = daeParameter("Dsld_{l}".format(l=l),
                    unit(), self,
                    "Diffusivity in electrode active particles",
                    [self.Nvol_ac[l], self.Npart_ac[l]])
            self.kappa_ac[l] = daeParameter("kappa_{l}".format(l=l),
                    unit(), self,
                    "kappa for each particle",
                    [self.Nvol_ac[l], self.Npart_ac[l]])
            self.Omga_ac[l] = daeParameter("Omga_{l}".format(l=l),
                    unit(), self,
                    "regular solution parameter for each particle [J]",
                    [self.Nvol_ac[l], self.Npart_ac[l]])
            self.Omgb_ac[l] = daeParameter("Omgb_{l}".format(l=l),
                    unit(), self,
                    "regular solution parameter for each particle [J]",
                    [self.Nvol_ac[l], self.Npart_ac[l]])
            self.Omgc_ac[l] = daeParameter("Omgc_{l}".format(l=l),
                    unit(), self,
                    "regular solution parameter for each particle [J]",
                    [self.Nvol_ac[l], self.Npart_ac[l]])
            self.k0_ac[l] = daeParameter("k0_{l}".format(l=l),
                    unit(), self,
                    "exchange current density rate constant for each particle",
                    [self.Nvol_ac[l], self.Npart_ac[l]])
            self.beta_s_ac[l] = daeParameter("beta_s_{l}".format(l=l),
                    unit(), self,
                    "surface wetting, nondim: kappa*d(gamma_s)/dc",
                    [self.Nvol_ac[l], self.Npart_ac[l]])
            self.delta_L_ac[l] = daeParameter("delta_L_{l}".format(l=l),
                    unit(), self,
                    "Length ratios for particle: Vp/(Ap*Rp)",
                    [self.Nvol_ac[l], self.Npart_ac[l]])
            self.MHC_Aa_ac[l] = daeParameter("MHC_Aa_{l}".format(l=l),
                    unit(), self,
                    "MHC factor for erf approximation parameter",
#                    [self.Ntrode])
                    [self.Nvol_ac[l], self.Npart_ac[l]])
            self.scond_ac[l] = daeParameter("scond_{l}".format(l=l),
                    unit(), self,
                    "surface conductivity of particles",
                    [self.Nvol_ac[l], self.Npart_ac[l]])
            self.G_ac[l] = daeParameter("G_{l}".format(l=l),
                    unit(), self,
                    "conductance particles",
                    [self.Nvol_ac[l], self.Npart_ac[l]])
            self.psd_num_ac[l] = daeParameter("psd_numVols_{l}".format(l=l),
                    unit(), self,
                    "Particle numbers of discretizations",
                    [self.Nvol_ac[l], self.Npart_ac[l]])
            self.psd_len_ac[l] = daeParameter("psd_lengths_{l}".format(l=l),
                    unit(), self,
                    "Particle lengths [nm]",
                    [self.Nvol_ac[l], self.Npart_ac[l]])
            self.psd_area_ac[l] = daeParameter("psd_active_areas_{l}".format(l=l),
                    unit(), self,
                    "Particle active areas [nm^2]",
                    [self.Nvol_ac[l], self.Npart_ac[l]])
            self.psd_vol_ac[l] = daeParameter("psd_volumes_{l}".format(l=l),
                    unit(), self,
                    "Particle volumes [nm^3]",
                    [self.Nvol_ac[l], self.Npart_ac[l]])
        self.dphi_eq_ref_ac = daeParameter("dphi_eq_ref_ac",
                unit(), self,
                "dimensionless potential offset in referencing fit " +
                "delta_phi_eq curves -- used only for initialization",
                [self.Ntrode])
        self.alpha_ac = daeParameter("alpha_ac", unit(), self,
                " Charge transfer coefficient",
                [self.Ntrode])
        self.dim_csmax_ac = daeParameter("dim_csmax_ac", unit(), self,
                "maximum lithium concentration in solid [mol/m^3]",
                [self.Ntrode])
        self.MHC_erfstretch_ac = daeParameter("MHC_b_ac", unit(), self,
                "MHC factor for erf approximation parameter",
                [self.Ntrode])
        self.T = daeParameter("T", unit(), self,
                "Non dimensional temperature")
        self.currset = daeParameter("currset", unit(), self,
                "dimensionless current")
        self.Vset = daeParameter("Vset", unit(), self,
                "dimensionless applied voltage (relative to " +
                "Delta V OCV of the  cell)")
        self.cwet_ac = daeParameter("cwet_ac", unit(), self,
                "Wetted surface concentration",
                [self.Ntrode])
        self.B_ac = daeParameter("B_ac", unit(), self,
                "Stress coefficient for each particle",
                [self.Ntrode])
        self.EvdW_ac = daeParameter("EvdW_ac", unit(), self,
                "vdW energy between planes",
                [self.Ntrode])
        self.lmbda_ac = daeParameter("lmbda_ac", unit(), self,
                "Marcus reorganizational energy",
                [self.Ntrode])
        self.mcond_ac = daeParameter("mcond_ac", unit(), self,
                "conductivity of cathode",
                [self.Ntrode])
        self.z = daeParameter("z", unit(), self,
                "capacity ratio of cathode:anode")

    def DeclareEquations(self):
        daeModel.DeclareEquations(self)

        # Some values of domain lengths
        trodes = self.trodes
        Nvol_ac = [0, 0]
        Npart_ac = [0, 0]
        for l in trodes:
            Nvol_ac[l] = self.Nvol_ac[l].NumberOfPoints
            Npart_ac[l] = self.Npart_ac[l].NumberOfPoints
        if self.D['Nvol_s'] >= 1: # if we have a separator
            Nvol_s = self.Nvol_s.NumberOfPoints
        else:
            Nvol_s = 0
        Nlyte = Nvol_s + np.sum(Nvol_ac)
        Nsld_mat_ac = np.empty(2, dtype=object)
        for l in trodes:
            Nsld_mat_ac[l] = np.zeros((Nvol_ac[l], Npart_ac[l]), dtype=np.integer)
            for i in range(Nvol_ac[l]):
                for j in range(Npart_ac[l]):
                    Nptsij = self.Nsld_mat_ac[l][i, j].NumberOfPoints
                    Nsld_mat_ac[l][i, j] = Nptsij
        # External function -- erf -- prepare to store external
        # function objects. For some reason, each external function
        # object that gets created has to stay 'alive' as an attribute
        # of the model, but we won't have to keep track of indexing.
        self.erfvec = []

#        # Prepare the noise
#        # maybe "numnoise" should be a parameter?
#        numnoise = 200/10.
#        noise_prefac = 2e-4
#        # a vector going from 0 to the max simulation time.
#        time_vec = np.linspace(0, (1./self.currset.GetValue()), numnoise)
#        # daeScalarExternalFunction (noise interpolation done as vector)
#        self.noise_local = np.empty(2, dtype=object)
#        for l in trodes:
#            self.noise_local[l] = np.empty((Nvol_ac[l], Npart_ac[l]),
#                    dtype=object)
#            for i in range(Nvol_ac[l]):
#                for j in range(Npart_ac[l]):
#                    Nsld = Nsld_mat_ac[l][i, j]
#                    noise_data = noise_prefac*np.random.randn(numnoise, Nsld)
#                    # Previous_output is common for all external functions
#                    previous_output = []
#                    self.noise_local[l][i, j] = [noise("Noise{l}{i}{j}".format(l=l,i=i,j=j), self,
#                        unit(), Time(), time_vec, noise_data,
#                        previous_output, _position_) for _position_
#                        in range(Nsld)]

        # Define the average concentration in each particle (algebraic
        # equations)
        for l in trodes:
            for i in range(Nvol_ac[l]):
                for j in range(Npart_ac[l]):
                    Nij = Nsld_mat_ac[l][i, j]
                    r_vec, volfrac_vec = self.get_unit_solid_discr(
                            self.D['solidShape_ac'][l],
                            self.D['solidType_ac'][l], Nij)
                    eq1 = self.CreateEquation("c1bar_trode{l}vol{i}part{j}".format(i=i,j=j,l=l))
                    eq2 = self.CreateEquation("c2bar_trode{l}vol{i}part{j}".format(i=i,j=j,l=l))
                    eq1.Residual = self.c1bar_sld_ac[l](i, j)
                    eq2.Residual = self.c2bar_sld_ac[l](i, j)
                    for k in range(Nij):
                        eq1.Residual -= self.c1_sld_ac[l][i, j](k) * volfrac_vec[k]
                        eq2.Residual -= self.c2_sld_ac[l][i, j](k) * volfrac_vec[k]
#                    eq.BuildJacobianExpressions = True

        # Define the overall filling fraction in the electrodes
        for l in trodes:
            eq = self.CreateEquation("ffrac_{l}".format(l=l))
            eq.Residual = self.ffrac_ac[l]()
            # Make a float of Vtot, total particle volume in electrode
            # Note: for some reason, even when "factored out", it's a bit
            # slower to use Sum(self.psd_vol_ac[l].array([], [])
            Vtot = np.sum(self.psd_vol_ac[l].npyValues)
            tmp = 0
            for i in range(Nvol_ac[l]):
                for j in range(Npart_ac[l]):
                    Vpart = self.psd_vol_ac[l](i, j)
                    # For some reason the following slower, so
                    # it's helpful to factor out the Vtot sum:
                    # eq.Residual -= self.cbar_sld_ac[l](i, j) * (Vpart/Vtot)
                    tmp += 0.5 * (self.c1bar_sld_ac[l](i, j) +
                            self.c2bar_sld_ac[l](i, j) ) * Vpart
            eq.Residual -= tmp / Vtot

        # Define dimensionless j_plus for each electrode volume
        for l in trodes:
            for i in range(Nvol_ac[l]):
                eq = self.CreateEquation("j_plus_trode{l}vol{i}".format(i=i,l=l))
                # Start with no reaction, then add reactions for each
                # particle in the volume.
                res = 0
                # sum over particle volumes in given electrode volume
                Vu = Sum(self.psd_vol_ac[l].array(i, []))
                for j in range(Npart_ac[l]):
                    # The volume of this particular particle
                    Vj = self.psd_vol_ac[l](i, j)
                    Nij = Nsld_mat_ac[l][i, j]
                    r_vec, volfrac_vec = self.get_unit_solid_discr(
                            self.D['solidShape_ac'][l],
                            self.D['solidType_ac'][l], Nij)
                    tmp = 0
                    for k in range(Nij):
                        tmp += (0.5 * (self.c1_sld_ac[l][i, j].dt(k) +
                                self.c2_sld_ac[l][i, j].dt(k)) *
                                volfrac_vec[k])
                    res += tmp * Vj
                eq.Residual = self.j_plus_ac[l](i) - res/Vu

        # Solid active particle concentrations, potential, and bulk
        # solid potential
        for l in trodes:
            # Solid active particles concentration/potential
            for i in range(Nvol_ac[l]):
                # Calculate the solid concentration rates of change
                for j in range(Npart_ac[l]):
                    # Prepare the RHS function
                    Nij = Nsld_mat_ac[l][i, j]
                    # Note Mmat is often the identity matrix...
                    (Mmat, RHS_c1_sld_ij, RHS_c2_sld_ij) = self.calc_sld_dcs_dt(l, i, j)
                    dc1dt_vec = np.empty(Nij, dtype=object)
                    dc2dt_vec = np.empty(Nij, dtype=object)
                    dc1dt_vec[0:Nij] = [self.c1_sld_ac[l][i, j].dt(k) for k in range(Nij)]
                    dc2dt_vec[0:Nij] = [self.c2_sld_ac[l][i, j].dt(k) for k in range(Nij)]
                    LHS1_vec = self.MX(Mmat, dc1dt_vec)
                    LHS2_vec = self.MX(Mmat, dc2dt_vec)
#                    noisevec = self.noise_local[l][i, j] # NOISE
                    # Set up equations: Mmat*dcdt = RHS
                    for k in range(Nij):
                        eq = self.CreateEquation(
                                "dc1sdt_trode{l}vol{i}part{j}discr{k}".format(
                                    i=i,j=j,k=k,l=l))
#                        eq.Residual = self.c_sld[i, j].dt(k) - RHS_c_sld_ij[k]
                        eq.Residual = LHS1_vec[k] - RHS_c1_sld_ij[k]
#                        eq.Residual = (LHS_vec[k] - RHS_c_sld_ij[k] +
#                                noisevec[k]()) # NOISE
                        eq = self.CreateEquation(
                                "dc2sdt_trode{l}vol{i}part{j}discr{k}".format(
                                    i=i,j=j,k=k,l=l))
                        eq.Residual = LHS2_vec[k] - RHS_c2_sld_ij[k]

                # Also calculate the potential drop along cathode
                # particle surfaces, if desired
                simSurfCathCond = self.D['simSurfCond_ac'][l]
                if simSurfCathCond:
                    # Conservation of charge in the solid particles with
                    # Ohm's Law
                    LHS = self.calc_part_surf_LHS(l, i, j)
                    k0_part = self.k0_ac[l](i, j)
                    for k in range(Nij):
                        eq = self.CreateEquation(
                                "charge_cons_trode{l}vol{i}part{j}discr{k}".format(
                                    i=i,j=j,k=k,l=l))
                        RHS = 0.5 * (self.c1_sld_ac[l][i, j].dt(k) +
                                self.c2_sld_ac[l][i, j].dt(k)) / k0_part
                        eq.Residual = LHS[k] - RHS

            # Simulate the potential drop along the bulk electrode
            # solid phase
            simBulkCond = self.D['simBulkCond_ac'][l]
            if simBulkCond:
                # Calculate the RHS for electrode conductivity
                phi_tmp = np.empty(Nvol_ac[l]+2, dtype=object)
                phi_tmp[1:-1] = [self.phi_ac[l](i) for i in range(Nvol_ac[l])]
                if l == 0: # anode
                    # Potential at the current collector is from
                    # simulation
                    phi_tmp[0] = self.phi_applied()
                    # No current passes into the electrolyte
                    phi_tmp[-1] = phi_tmp[-2]
                else: # cathode
                    phi_tmp[0] = phi_tmp[1]
                    # Potential at current at current collector is
                    # reference (set)
                    phi_tmp[-1] = self.phi_cathode()
                dx = 1./Nvol_ac[l]
                RHS_phi_tmp = -np.diff(-self.mcond_ac(l)*np.diff(phi_tmp)/dx)/dx
            # Actually set up the equations for bulk solid phi
            for i in range(Nvol_ac[l]):
                eq = self.CreateEquation("phi_ac_trode{l}vol{i}".format(i=i,l=l))
                if simBulkCond:
                    eq.Residual = (-self.epsbeta_ac(l)*self.j_plus_ac[l](i) -
                            RHS_phi_tmp[i])
                else:
                    if l == 0: # anode
                        eq.Residual = self.phi_ac[l](i) - self.phi_applied()
                    else: # cathode
                        eq.Residual = self.phi_ac[l](i) - self.phi_cathode()

            # Simulate the potential drop along the connected
            # particles
            simPartCond = self.D['simPartCond_ac'][l]
            for i in range(Nvol_ac[l]):
                phi_bulk = self.phi_ac[l](i)
                for j in range(Npart_ac[l]):
                    G_l = self.G_ac[l](i, j)
                    phi_n = self.phipart_ac[l](i, j)
                    if j == 0: # reference bulk phi
                        phi_l = self.phi_ac[l](i)
                    else:
                        phi_l = self.phipart_ac[l](i, j - 1)
                    if j == (Npart_ac[l] - 1): # No particle at and of "chain"
                        G_r = 0
                        phi_r = phi_n
                    else:
                        G_r = self.G_ac[l](i, j + 1)
                        phi_r = self.phipart_ac[l](i, j + 1)
                    # Find average dcs/dt for this particle
                    Nij = Nsld_mat_ac[l][i, j]
                    r_vec, volfrac_vec = self.get_unit_solid_discr(
                            self.D['solidShape_ac'][l],
                            self.D['solidType_ac'][l], Nij)
                    dcsbardt = 0
                    for k in range(Nij):
                        dcsbardt += (self.c1_sld_ac[l][i, j].dt(k) *
                                volfrac_vec[k])
                        dcsbardt += (self.c2_sld_ac[l][i, j].dt(k) *
                                volfrac_vec[k])
                    # charge conservation equation around this particle
                    eq = self.CreateEquation("phi_ac_trode{l}vol{i}part{j}".format(i=i,l=l,j=j))
                    if simPartCond:
                        # -dcsbar/dt = I_l - I_r
                        eq.Residual = dcsbardt + (
                                (-G_l * (phi_n - phi_l)) -
                                (-G_r * (phi_r - phi_n)))
                    else:
                        eq.Residual = self.phipart_ac[l](i, j) - self.phi_ac[l](i)

        # If we have a single electrode volume (in a perfect bath),
        # electrolyte equations are simple
        if Nvol_ac[0] == 0 and Nvol_s == 0 and Nvol_ac[1] == 1:
            eq = self.CreateEquation("c_lyte")
            eq.Residual = self.c_lyte_ac[1].dt(0) - 0
            eq = self.CreateEquation("phi_lyte")
            eq.Residual = self.phi_lyte_ac[1](0) - self.phi_applied()
        else:
            # Calculate RHS for electrolyte equations
            c_lyte = np.empty(Nlyte, dtype=object)
            c_lyte[0:Nvol_ac[0]] = [self.c_lyte_ac[0](i)
                    for i in range(Nvol_ac[0])] # anode
            c_lyte[Nvol_ac[0]:Nvol_ac[0] + Nvol_s] = [self.c_lyte_s(i)
                    for i in range(Nvol_s)] # separator
            c_lyte[Nvol_ac[0] + Nvol_s:Nlyte] = [self.c_lyte_ac[1](i) for i in
                    range(Nvol_ac[1])] # cathode
            phi_lyte = np.empty(Nlyte, dtype=object)
            phi_lyte[0:Nvol_ac[0]] = [self.phi_lyte_ac[0](i)
                    for i in range(Nvol_ac[0])] # anode
            phi_lyte[Nvol_ac[0]:Nvol_ac[0] + Nvol_s] = [self.phi_lyte_s(i)
                    for i in range(Nvol_s)] # separator
            phi_lyte[Nvol_ac[0] + Nvol_s:Nlyte] = [self.phi_lyte_ac[1](i)
                    for i in range(Nvol_ac[1])] # cathode
            (RHS_c, RHS_phi) = self.calc_lyte_RHS(c_lyte, phi_lyte,
                    Nvol_ac, Nvol_s, Nlyte)
            # Equations governing the electrolyte in the separator
            offset = Nvol_ac[0]
            for i in range(Nvol_s):
                # Mass Conservation
                eq = self.CreateEquation(
                        "sep_lyte_mass_cons_vol{i}".format(i=i))
                eq.Residual = (self.poros_s()*self.c_lyte_s.dt(i) -
                        RHS_c[offset + i])
                # Charge Conservation
                eq = self.CreateEquation(
                        "sep_lyte_charge_cons_vol{i}".format(i=i))
                eq.Residual = (RHS_phi[offset + i])
            # Equations governing the electrolyte in the electrodes.
            # Here, we are coupled to the total reaction rates in the
            # solids.
            for l in trodes:
                if l == 0: # anode
                    offset = 0
                else: # cathode
                    offset = Nvol_ac[0] + Nvol_s
                for i in range(Nvol_ac[l]):
                    # Mass Conservation
                    eq = self.CreateEquation(
                            "lyteMassCons_trode{l}vol{i}".format(i=i,l=l))
                    eq.Residual = (self.poros_ac(l)*self.c_lyte_ac[l].dt(i) +
                            self.epsbeta_ac(l)*(1-self.tp())*self.j_plus_ac[l](i) -
                            RHS_c[offset + i])
                    # Charge Conservation
                    eq = self.CreateEquation(
                            "lyteChargeCons_trode{l}vol{i}".format(i=i,l=l))
                    eq.Residual = (self.epsbeta_ac(l)*self.j_plus_ac[l](i) -
                            RHS_phi[offset + i])

        # Define the total current. This can be done in either anode
        # or cathode equivalently.
        eq = self.CreateEquation("Total_Current")
        eq.Residual = self.current()
        limtrode = (1 if self.z.GetValue() < 1 else 0)
        dx = 1./Nvol_ac[limtrode]
        for i in range(Nvol_ac[limtrode]):
            # Factor of -1 is to get sign right if referencing anode
            eq.Residual -= dx * (-1)**(1-limtrode) * self.j_plus_ac[limtrode](i)

        Dp = self.D['Dp']
        Dm = self.D['Dm']
        zp = self.D['zp']
        zm = self.D['zm']
        Damb = ((zp+zm)*Dp*Dm)/(zp*Dp+zm*Dm)
        td = self.D['L_ac'][1]**2 / Damb
        if self.profileType == "CC":
            # Total Current Constraint Equation
            eq = self.CreateEquation("Total_Current_Constraint")
            if self.currset.GetValue() != 0.0:
                timeHorizon = 1./Abs(self.currset())
            else:
                timeHorizon = self.D['tend']/td
            eq.Residual = self.current() - self.currset()*(1 -
                    np.exp(-Time()/(timeHorizon*1e-3)))
        elif self.profileType == "CV":
            # Keep applied potential constant
            eq = self.CreateEquation("applied_potential")
            timeHorizon = self.D['tend']/td
            eq.Residual = self.phi_applied() - self.Vset()*(
#                    1)
#                    1 - np.exp(-Time()/(timeHorizon*1e-3)))
                    np.tanh(Time()/(45.0)))

        for eq in self.Equations:
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

    def calc_sld_dcs_dt(self, trode_indx, vol_indx, part_indx):
        # shorthand
        l = trode_indx
        i = vol_indx
        j = part_indx
        # Get some useful information
        simSurfCond = self.D['simSurfCond_ac'][l]
        solidType = self.D['solidType_ac'][l]
        solidShape = self.D['solidShape_ac'][l]
        rxnType = self.D['rxnType_ac'][l]
        delPhiEqFit = self.D['delPhiEqFit_ac'][l]
        # Get variables for this particle/electrode volume
        phi_lyte = self.phi_lyte_ac[l](i)
        phi_m = self.phipart_ac[l](i, j)
        c_lyte = self.c_lyte_ac[l](i)
        # Get the relevant parameters for this particle
        k0 = self.k0_ac[l](i, j)
        kappa = self.kappa_ac[l](i, j) # only used for ACR/CHR
        B = self.B_ac(l) # only used for ACR/CHR
        c1bar = self.c1bar_sld_ac[l](i, j) # only used for ACR/CHR
        c2bar = self.c2bar_sld_ac[l](i, j) # only used for ACR/CHR
        lmbda = self.lmbda_ac(l) # Only used for Marcus/MHC
        Aa = self.MHC_Aa_ac[l](i, j) # Only used for MHC
        b = self.MHC_erfstretch_ac(l) # Only used for MHC
        alpha = self.alpha_ac(l) # Only used for BV
        Omga = self.Omga_ac[l](i, j)
        Omgb = self.Omgb_ac[l](i, j)
        Omgc = self.Omgc_ac[l](i, j)
        Ds = self.Dsld_ac[l](i, j) # Only used for "diffn"
        # We need the (non-dimensional) temperature to get the
        # reaction rate dependence correct
        T = self.T()
        # Number of volumes in current particle
        Nij = self.Nsld_mat_ac[l][i, j].NumberOfPoints
        # Concentration (profile?) in the solid
        c1_sld = np.empty(Nij, dtype=object)
        c1_sld[:] = [self.c1_sld_ac[l][i, j](k) for k in range(Nij)]
        c2_sld = np.empty(Nij, dtype=object)
        c2_sld[:] = [self.c2_sld_ac[l][i, j](k) for k in range(Nij)]
        # Assume dilute electrolyte
        act_O = c_lyte
        mu_O = T*np.log(Max(eps, act_O))
        # Calculate chemical potential of reduced state

        if solidType in ["ACR", "homog", "homog_sdn"]:
            if solidType == "ACR":
                # Make a blank array to allow for boundary conditions
                mu_R = self.calc_ACR_mu_R(c_sld, Omga, B, kappa, T)
                # If we're also simulating potential drop along the solid,
                # use that instead of self.phi_c(i)
                if simSurfCond:
                    phi_m = np.empty(Nij, dtype=object)
                    phi_m[:] = [self.phi_sld_ac[l][i, j](k) for k in range(Nij)]
            elif solidType in ["homog", "homog_sdn"]:
                mu_R = self.mu_reg_sln(c_sld, Omga, T, eps)
            # XXX -- Temp dependence!
            act_R = np.exp(mu_R/T)
            # eta = electrochem pot_R - electrochem pot_O
            # eta = (mu_R + phi_R) - (mu_O + phi_O)
            if delPhiEqFit:
                material = self.D['material_ac'][l]
                fits = delta_phi_fits.DPhiFits(self.D)
                phifunc = fits.materialData[material]
                delta_phi_eq = phifunc(c_sld[-1], self.dphi_eq_ref_ac(l))
                eta = (phi_m - phi_lyte) - delta_phi_eq
            else:
                eta = (mu_R + phi_m) - (mu_O + phi_lyte)
            if rxnType == "Marcus":
                Rate = self.R_Marcus(k0, lmbda, c_lyte, c_sld, eta, T)
            elif rxnType == "BV":
                Rate = self.R_BV(k0, alpha, c_sld, act_O, act_R, eta, T)
            elif rxnType == "MHC":
#                Rate = self.R_MHC(k0, lmbda, eta, Aa, b, T)
                Rate = self.R_MHC(k0, lmbda, eta, Aa, b, T, c_sld)
            M = sprs.eye(Nij, Nij, format="csr")
            return (M, Rate)

        elif solidType in ["diffn", "CHR"] and solidShape in ["sphere", "cylinder"]:
            # For discretization background, see Zeng & Bazant 2013
            # Mass matrix is common for spherical shape, diffn or CHR
            Rs = 1. # (non-dimensionalized by itself)
            r_vec, volfrac_vec = self.get_unit_solid_discr(
                    solidShape, solidType, Nij)
            edges = np.hstack((0, (r_vec[0:-1] + r_vec[1:])/2, Rs))
            if solidShape == "sphere":
                Vp = 4./3. * np.pi * Rs**3
            elif solidShape == "cylinder":
                Vp = np.pi * Rs**2  # per unit height
            vol_vec = Vp * volfrac_vec
            dr = r_vec[1] - r_vec[0]
            M1 = sprs.diags([1./8, 3./4, 1./8], [-1, 0, 1],
                    shape=(Nij,Nij), format="csr")
            M1[1, 0] = M1[-2, -1] = 1./4
            M2 = sprs.diags(vol_vec, 0, format="csr")
            if solidShape == "sphere":
                M = M1*M2
            elif solidShape == "cylinder":
                M = M2

            # Diffn
            if solidType in ["diffn"]:
                if solidShape in ["cylinder"]:
                    raise NotImplementedError("cylinder-diffn!")
                # Take the surface concentration
                c_surf = c_sld[-1]
                # Assuming dilute solid
                act_R_surf = c_surf
                mu_R_surf = T*np.log(Max(eps, act_R_surf))
            # CHR
            elif solidType in ["CHR"]:
                # mu_R is like for ACR, except kappa*curv term
                EvdW = self.EvdW_ac(l)
                beta_s = self.beta_s_ac[l](i, j)
                mu1_R, mu2_R = self.calc_mu12_R(c1_sld, c2_sld, c1bar,
                        c2bar, Omga, Omgb, Omgc, B, kappa, EvdW, beta_s, T)
                act1_R = np.exp(mu1_R/T)
                act2_R = np.exp(mu2_R/T)
                mu1_R_surf = mu1_R[-1]
                mu2_R_surf = mu2_R[-1]
                # Take the surface concentration
                c1_surf = c1_sld[-1]
                c2_surf = c2_sld[-1]
                act1_R_surf = np.exp(mu1_R_surf/T)
                act2_R_surf = np.exp(mu2_R_surf/T)
            # Figure out overpotential
            if delPhiEqFit:
                material = self.D['material_ac'][l]
                fits = delta_phi_fits.DPhiFits(self.D)
                phifunc = fits.materialData[material]
                delta_phi_eq1 = phifunc(c1_surf, self.dphi_eq_ref_ac(l))
                delta_phi_eq2 = phifunc(c2_surf, self.dphi_eq_ref_ac(l))
                eta1 = (phi_m - phi_lyte) - delta_phi_eq1
                eta2 = (phi_m - phi_lyte) - delta_phi_eq2
            else:
                eta1 = (mu1_R_surf + phi_m) - (mu_O + phi_lyte)
                eta2 = (mu2_R_surf + phi_m) - (mu_O + phi_lyte)
            # Calculate reaction rate
            if rxnType == "Marcus":
                Rxn1 = self.R_Marcus(k0, lmbda, c_lyte, c1_surf, eta1, T)
                Rxn2 = self.R_Marcus(k0, lmbda, c_lyte, c2_surf, eta2, T)
            elif rxnType == "BV":
                Rxn1 = self.R_BV(k0, alpha, c1_surf, c2_surf, Omgb, act_O, act1_R_surf, eta1, T)
                Rxn2 = self.R_BV(k0, alpha, c2_surf, c1_surf, Omgb, act_O, act2_R_surf, eta2, T)
#                Rxn1 = self.R_BV(k0, alpha, c1_surf, act_O, act1_R_surf, eta1, T)
#                Rxn2 = self.R_BV(k0, alpha, c2_surf, act_O, act2_R_surf, eta2, T)
            elif rxnType == "MHC":
                Rxn1 = self.R_MHC(k0, lmbda, eta1, Aa, b, T, c1_surf)
                Rxn2 = self.R_MHC(k0, lmbda, eta2, Aa, b, T, c2_surf)
            # Finish up RHS discretization at particle surface
            if solidType in ["diffn"]:
#                RHS = np.empty(Nij, dtype=object)
                c_diffs = np.diff(c_sld)
#                RHS[1:Nij - 1] = 4*np.pi*(
#                        Ds*(r_vec[1:Nij - 1] + dr/2.)**2*c_diffs[1:]/dr -
#                        Ds*(r_vec[1:Nij - 1] - dr/2.)**2*c_diffs[:-1]/dr )
#                RHS[0] = 4*np.pi*Ds*(dr/2.)**2*c_diffs[0]/dr
                Flux_vec = self.calc_Flux_diffn(c_sld, Ds, Rxn, r_vec, dr)
                RHS = 4*np.pi*Flux_vec
#                RHS[-1] = 4*np.pi*(Rs**2 * self.delta_L_ac[l](i, j) * Rxn -
#                        Ds*(Rs - dr/2)**2*c_diffs[-1]/dr )
            elif solidType in ["CHR"]:
                Flux1_bc = 0.5*self.delta_L_ac[l](i, j)*Rxn1
                Flux2_bc = 0.5*self.delta_L_ac[l](i, j)*Rxn2
                # With chem potentials, can now calculate fluxes
                Flux1_vec, Flux2_vec = self.calc_Flux12(c1_sld,
                        c2_sld, mu1_R, mu2_R, Ds, Flux1_bc, Flux2_bc, dr, T)
                if solidShape == "sphere":
                    area_vec = 4*np.pi*edges**2
                elif solidShape == "cylinder":
                    area_vec = 2*np.pi*edges  # per unit height
                FluxTerm1 = np.diff(Flux1_vec*area_vec)
                FluxTerm2 = np.diff(Flux2_vec*area_vec)
#                kinterlayer = 1e2
#                interLayerRxn = (kinterlayer * (1 - c1_sld) *
#                        (1 - c2_sld) * (act1_R - act2_R))
#                RxnTerm1 = -interLayerRxn
#                RxnTerm2 = interLayerRxn
                RxnTerm1 = 0
                RxnTerm2 = 0
                RHS1 = FluxTerm1 + RxnTerm1
                RHS2 = FluxTerm2 + RxnTerm2
            return (M, RHS1, RHS2)

    def calc_ACR_mu_R(self, c_sld, Omga, B, kappa, T, l):
        Nij = len(c_sld)
        cstmp = np.empty(Nij+2, dtype=object)
        cstmp[1:-1] = c_sld
        cstmp[0] = self.cwet_ac(l)
        cstmp[-1] = self.cwet_ac(l)
        dxs = 1./Nij
        curv = np.diff(cstmp, 2)/(dxs**2)
        mu_R = ( self.mu_reg_sln(c_sld, Omga, T, eps) - kappa*curv
                + self.B_ac(l)*(c_sld - cbar) )
        return mu_R

    def calc_mu12_R(self, c1_sld, c2_sld, c1bar, c2bar,
            Omga, Omgb, Omgc, B, kappa, EvdW, beta_s, T):
        Nij = len(c1_sld)
        solidType = "CHR"
        solidShape = "cylinder"
        Rs = 1. # (non-dimensionalized by itself)
        r_vec, volfrac_vec = self.get_unit_solid_discr(
                solidShape, solidType, Nij)
        dr = r_vec[1] - r_vec[0]
        # mu_R is like for ACR, except kappa*curv term
        mu1_R = ( self.mu_reg_sln(c1_sld, Omga, T, eps) +
                B*(c1_sld - c1bar) )
        mu2_R = ( self.mu_reg_sln(c2_sld, Omga, T, eps) +
                B*(c2_sld - c2bar) )
        mu1_reg = self.mu_reg_sln(c1_sld, Omga, T, eps)
        mu2_reg = self.mu_reg_sln(c2_sld, Omga, T, eps)
        mu1_B = B*(c1_sld - c1bar)
        mu2_B = B*(c2_sld - c2bar)
        mu1_vdW = EvdW*(30*c1_sld**2*(1-c1_sld)**2)
        mu2_vdW = EvdW*(30*c2_sld**2*(1-c2_sld)**2)
        mu1_intr = Omgb*c2_sld + Omgc*c2_sld*(1-c2_sld)*(1-2*c1_sld)
        mu2_intr = Omgb*c1_sld + Omgc*c1_sld*(1-c1_sld)*(1-2*c2_sld)
        # Graphite vdW interactions?
#        EvdW = self.EvdW_ac(l)
        mu1_R += EvdW*(30*c1_sld**2*(1-c1_sld)**2)
        mu2_R += EvdW*(30*c2_sld**2*(1-c2_sld)**2)
        # Surface conc gradient given by natural BC
#        beta_s = self.beta_s_ac[l](i, j)
        if solidShape == "sphere":
            raise NotImplementedError("no sphere in this version")
            mu_R[0] -= 3 * kappa * (2*c_sld[1] - 2*c_sld[0])/dr**2
            mu_R[1:Nij - 1] -= kappa * (np.diff(c_sld, 2)/dr**2 +
                    (c_sld[2:] - c_sld[0:-2])/(dr*r_vec[1:-1])
                    )
            mu_R[Nij - 1] -= kappa * ((2./Rs)*beta_s +
                    (2*c_sld[-2] - 2*c_sld[-1] + 2*dr*beta_s)/dr**2
                    )
        elif solidShape == "cylinder":
            mu1_R[0] -= 2 * kappa * (2*c1_sld[1] - 2*c1_sld[0])/dr**2
            mu2_R[0] -= 2 * kappa * (2*c2_sld[1] - 2*c2_sld[0])/dr**2
            mu1_R[1:Nij - 1] -= kappa * (np.diff(c1_sld, 2)/dr**2 +
                    (c1_sld[2:] - c1_sld[0:-2])/(2 * dr*r_vec[1:-1])
                    )
            mu2_R[1:Nij - 1] -= kappa * (np.diff(c2_sld, 2)/dr**2 +
                    (c2_sld[2:] - c2_sld[0:-2])/(2 * dr*r_vec[1:-1])
                    )
            mu1_R[Nij - 1] -= kappa * ((1./Rs)*beta_s +
                    (2*c1_sld[-2] - 2*c1_sld[-1] + 2*dr*beta_s)/dr**2
                    )
            mu2_R[Nij - 1] -= kappa * ((1./Rs)*beta_s +
                    (2*c2_sld[-2] - 2*c2_sld[-1] + 2*dr*beta_s)/dr**2
                    )
            mu1_R += Omgb*c2_sld + Omgc*c2_sld*(1-c2_sld)*(1-2*c1_sld)
            mu2_R += Omgb*c1_sld + Omgc*c1_sld*(1-c1_sld)*(1-2*c2_sld)
        if (type(c1_sld[0]) == pyCore.adouble):
            return (mu1_R, mu2_R)
        else:
            return (mu1_R, mu2_R, mu1_reg, mu2_reg, mu1_B, mu2_B, mu1_vdW,
                mu2_vdW, mu1_intr, mu2_intr)
#        return (mu1_R, mu2_R)
#        return (mu1_R, mu2_R, mu1_reg, mu2_reg, mu1_B, mu2_B, mu1_vdW,
#                mu2_vdW, mu1_intr, mu2_intr)

    def calc_Flux_diffn(self, c_sld, Ds, Rxn, r_vec, dr):
        Nij = len(c_sld)
        c_diffs = np.diff(c_sld)
        Flux_vec = np.empty(Nij, dtype=object)
        Flux_vec[1:Nij-1] =  ( Ds*(r_vec[1:Nij - 1] + dr/2.)**2*c_diffs[1:]/dr -
                Ds*(r_vec[1:Nij - 1] - dr/2.)**2*c_diffs[:-1]/dr )
        Flux_vec[0] = Ds*(dr/2.)**2*c_diffs[0]/dr
        Flux_vec[-1] = ( Rs**2 * self.delta_L_ac[l](i, j) * Rxn -
                        Ds*(Rs - dr/2)**2*c_diffs[-1]/dr )
        return Flux_vec

    def calc_Flux12(self, c1_sld, c2_sld, mu1_R, mu2_R, Ds, Flux1_bc,
            Flux2_bc, dr, T):
        Nij = len(c1_sld)
        Flux1_vec = np.empty(Nij+1, dtype=object)
        Flux2_vec = np.empty(Nij+1, dtype=object)
        Flux1_vec[0] = 0  # Symmetry at r=0
        Flux2_vec[0] = 0  # Symmetry at r=0
        Flux1_vec[-1] = Flux1_bc
        Flux2_vec[-1] = Flux2_bc
        c1_edges = (c1_sld[0:-1]  + c1_sld[1:])/2.
        c2_edges = (c2_sld[0:-1]  + c2_sld[1:])/2.
        cbar_edges = 0.5*(c1_edges + c2_edges)
        Flux1_vec[1:Nij] = (Ds/T * (1 - c1_edges)**(1.00) * c1_edges *
                np.diff(mu1_R)/dr)
        Flux2_vec[1:Nij] = (Ds/T * (1 - c2_edges)**(1.00) * c2_edges *
                np.diff(mu2_R)/dr)
#        mu1_edges = (mu1_R[0:-1]  + mu1_R[1:])/2.
#        mu2_edges = (mu2_R[0:-1]  + mu2_R[1:])/2.
#        gamma1_TS_edges = np.exp(mu1_edges/T) / c1_edges
#        gamma2_TS_edges = np.exp(mu2_edges/T) / c2_edges
#        Flux1_vec[1:Nij] = (Ds/T * c1_edges / gamma1_TS_edges *
#                np.diff(mu1_R)/dr)
#        Flux2_vec[1:Nij] = (Ds/T * c2_edges / gamma2_TS_edges *
#                np.diff(mu2_R)/dr)
        return (Flux1_vec, Flux2_vec)

    def calc_lyte_RHS(self, cvec, phivec, Nvol_ac, Nvol_s, Nlyte):
        # Discretization
        # The lengths are nondimensionalized by the cathode length
        dxvec = np.empty(np.sum(Nvol_ac) + Nvol_s + 2, dtype=object)
        dxa = np.array(Nvol_ac[0] * [self.L_ac(0)/Nvol_ac[0]])
        dxc = np.array(Nvol_ac[1] * [self.L_ac(1)/Nvol_ac[1]])
        dxs = np.array(Nvol_s * [self.L_s()/Nvol_s])
        dxtmp = np.hstack((dxa, dxs, dxc))
        dxvec[1:-1] = dxtmp
        dxvec[0] = dxvec[1]
        dxvec[-1] = dxvec[-2]
        dxd1 = (dxvec[0:-1] + dxvec[1:]) / 2.
        dxd2 = dxtmp

        # The porosity vector
        porosvec = np.empty(Nlyte + 1, dtype=object)
        # Use the Bruggeman relationship to approximate an effective
        # effect on the transport.
        porosvec[0:Nvol_ac[0]] = [self.poros_ac(0)**(3./2)
            for i in range(Nvol_ac[0])] # anode
        porosvec[Nvol_ac[0]:Nvol_ac[0] + Nvol_s] = [self.poros_s()**(3./2)
            for i in range(Nvol_s)] # separator
        porosvec[Nvol_ac[0] + Nvol_s:Nlyte+1] = [self.poros_ac(1)**(3./2)
            for i in range(Nvol_ac[1]+1)] # cathode

        # Mass conservation equations
        ctmp = np.empty(Nlyte + 2, dtype=object)
        ctmp[1:-1] = cvec
        # If we don't have a real anode, the total current flowing
        # into the electrolyte is set
        limtrode = (1 if self.z.GetValue() < 1 else 0)
        if Nvol_ac[0] == 0:
            ctmp[0] = ( ctmp[1] + (self.current() *
                self.epsbeta_ac(limtrode) *
                (1-self.tp())*dxvec[0])/porosvec[0] )
        else: # porous anode -- no elyte flux at anode current collector
            ctmp[0] = ctmp[1]
        # No electrolyte flux at the cathode current collector
        ctmp[-1] = ctmp[-2]
        # Diffusive flux in the electrolyte
        cflux = -porosvec*np.diff(ctmp)/dxd1
        # Divergence of the flux
        RHS_c = -np.diff(cflux)/dxd2

        # Charge conservation equations
        phitmp = np.empty(Nlyte + 2, dtype=object)
        phitmp[1:-1] = phivec
        # If we don't have a full anode, assume no rxn resistance at a
        # lithium anode, and measure relative to Li
        if Nvol_ac[0] == 0:
            phitmp[0] = self.phi_applied()
        else: # porous anode -- no flux into anode current collector
            phitmp[0] = phitmp[1]
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
        currdens = (-((Dp - Dm)*np.diff(ctmp)/dxd1) -
                (zp*Dp + zm*Dm)*c_edges*np.diff(phitmp)/dxd1)
        RHS_phi = -np.diff(porosvec*currdens)/dxd2
        return (RHS_c, RHS_phi)

    def calc_part_surf_LHS(self, trode_indx, vol_indx, part_indx):
        # shorthand
        l = trode_indx
        i = vol_indx
        j = part_indx
        # Number of volumes in current particle
        Nij = self.Nsld_mat_ac[l][i, j].NumberOfPoints
        # solid potential variables for this particle
        phi_tmp = np.empty(Nij + 2, dtype=object)
        phi_tmp[1:-1] = [self.phi_sld_ac[l][i, j](k) for k in
                range(Nij)]
        # BC's -- "touching carbon at each end"
        phi_s_local = self.phipart_ac[l](i, j)
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

    def mu_reg_sln(self, c, Omga, T, eps):
        if (type(c[0]) == pyCore.adouble):
            maxfunc = Max
        else:
            maxfunc = max
        return np.array([ Omga*(1-2*c[i])
                + T*np.log(maxfunc(eps, c[i])/maxfunc(eps, 1-c[i]))
                for i in range(len(c)) ])

    def R_BV(self, k0, alpha, ci_sld, cj_sld, Omgb, act_O, act_R, eta, T):
        gamma_ts = (1./(ci_sld**1.0*(1-ci_sld)))#/ci_sld**(0.5)#*np.exp((Omgb/T)*cj_sld)
#    def R_BV(self, k0, alpha, c_sld, act_O, act_R, eta, T):
#        gamma_ts = (1./(1-c_sld))
        ecd = ( k0 * act_O**(1-alpha)
                * act_R**(alpha) / gamma_ts )
        Rate = ( ecd *
            (np.exp(-alpha*eta/T) - np.exp((1-alpha)*eta/T)) )
        return Rate

    def R_Marcus(self, k0, lmbda, c_lyte, c_sld, eta, T):
        c_sld = np.array([Max(eps, c_sld[i]) for i in
            range(len(c_sld))])
        alpha = 0.5*(1 + (T/lmbda) * np.log(Max(eps, c_lyte)/c_sld))
        # We'll assume c_e = 1 (at the standard state for electrons)
#        ecd = ( k0 * np.exp(-lmbda/(4.*T)) *
#        ecd = ( k0 *
        ecd = ( k0 * (1-c_sld) *
                c_lyte**((3-2*alpha)/4.) *
                c_sld**((1+2*alpha)/4.) )
        Rate = ( ecd * np.exp(-eta**2/(4.*T*lmbda)) *
            (np.exp(-alpha*eta/T) - np.exp((1-alpha)*eta/T)) )
        return Rate

#    def R_MHC(self, k0, lmbda, eta, Aa, b, T):
    def R_MHC(self, k0, lmbda, eta, Aa, b, T, c_sld):
        def knet(eta, lmbda, b):
            argox = (eta - lmbda)/b
            argrd = (-eta - lmbda)/b
            k = Aa*((Erf(argrd) + 1) - (Erf(argox) +1))
            return k
        if type(eta) == np.ndarray:
            Rate = np.empty(len(eta), dtype=object)
            for i, etaval in enumerate(eta):
#                Rate[i] = knet(etaval, lmbda, b)
                Rate[i] = (1-c_sld[i])*knet(etaval, lmbda, b)
        else:
#            Rate = np.array([knet(eta, lmbda, b)])
            Rate = np.array([(1-c_sld)*knet(eta, lmbda, b)])
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

    def get_unit_solid_discr(self, solidShape, solidType, Nij):
        if solidShape == "C3" and solidType in ["ACR"]:
            r_vec = None
            # For 1D particle, the vol fracs are simply related to the
            # length discretization
            volfrac_vec = (1./Nij) * np.ones(Nij)  # scaled to 1D particle volume
            return r_vec, volfrac_vec
        if solidType in ["homog", "homog_sdn"]:
            r_vec = None
            volfrac_vec = np.ones(1)
            return r_vec, volfrac_vec
        if solidShape == "sphere" and solidType in ["diffn", "CHR"]:
            Rs = 1.
            dr = Rs/(Nij - 1)
            r_vec = np.linspace(0, Rs, Nij)
            vol_vec = 4*np.pi*(r_vec**2 * dr + (1./12)*dr**3)
            vol_vec[0] = 4*np.pi*(1./24)*dr**3
            vol_vec[-1] = (4./3)*np.pi*(Rs**3 - (Rs - dr/2.)**3)
            Vp = 4./3.*np.pi*Rs**3
            volfrac_vec = vol_vec/Vp
            return r_vec, volfrac_vec
        if solidShape == "cylinder" and solidType in ["diffn", "CHR"]:
            Rs = 1.
            h = 1.
            dr = Rs / (Nij - 1)
            r_vec = np.linspace(0, Rs, Nij)
            vol_vec = np.pi * h * 2 * r_vec * dr
            vol_vec[0] = np.pi * h * dr**2 / 4.
            vol_vec[-1] = np.pi * h * (Rs * dr - dr**2 / 4.)
            Vp = np.pi * Rs**2 * h
            volfrac_vec = vol_vec / Vp
            return r_vec, volfrac_vec
        else:
            raise NotImplementedError("Fix shape volumes!")

class simMPET(daeSimulation):
    def __init__(self, D=None):
        daeSimulation.__init__(self)
        if D is None:
            raise Exception("Need parameters input")
        self.D = D
        # TODO -- why is psd generation in __init__??
        self.Nvol_s = D['Nvol_s']
        self.Nvol_ac = [D['Nvol_ac'][0], D['Nvol_ac'][1]]
        self.Npart_ac = [D['Npart_ac'][0], D['Npart_ac'][1]]
        if self.Nvol_ac[0] >= 1: # If we have a full anode
            self.trodes = [0, 1]
        else:
            self.trodes = [1]
        self.test_input(D)

        # Generate psd, conductivity distribution
        self.psd_num = np.empty(2, dtype=object)
        self.psd_len = np.empty(2, dtype=object)
        self.G = np.empty(2, dtype=object)
        for l in self.trodes:
            psd_mean = D['psd_mean_ac'][l]
            psd_stddev = D['psd_stddev_ac'][l]
            G_mean = D['G_mean_ac'][l]
            G_stddev = D['G_stddev_ac'][l]
            Nvol = self.Nvol_ac[l]
            Npart = self.Npart_ac[l]
            solidType = D['solidType_ac'][l]
            # Make a length-sampled particle size distribution
            # Log-normally distributed
            if psd_stddev == 0:
                psd_raw = psd_mean*np.ones((Nvol, Npart))
            else:
                psd_var = psd_stddev**2
                mu = np.log((psd_mean**2)/np.sqrt(psd_var+psd_mean**2))
                sigma = np.sqrt(np.log(psd_var/(psd_mean**2)+1))
                psd_raw = np.random.lognormal(mu, sigma,
                        size=(Nvol, Npart))
            # For particles with internal profiles, convert psd to
            # integers -- number of steps
            if solidType in ["ACR"]:
                solidDisc = D['solidDisc_ac'][l]
                self.psd_num[l] = np.ceil(psd_raw/solidDisc).astype(np.integer)
                self.psd_len[l] = solidDisc*self.psd_num[l]
            elif solidType in ["CHR", "diffn"]:
                solidDisc = D['solidDisc_ac'][l]
                self.psd_num[l] = np.ceil(psd_raw/solidDisc).astype(np.integer) + 1
                self.psd_len[l] = solidDisc*(self.psd_num[l] - 1)
            # For homogeneous particles (only one "volume" per particle)
            elif solidType in ["homog", "homog_sdn"]:
                # Each particle is only one volume
                self.psd_num[l] = np.ones(psd_raw.shape).astype(np.integer)
                # The lengths are given by the original length distr.
                self.psd_len[l] = psd_raw
            else:
                raise NotImplementedError("Solid types missing here")

            if G_stddev == 0:
                self.G[l] = G_mean * np.ones((Nvol, Npart))
            else:
                G_var = G_stddev**2
                mu = np.log((G_mean**2)/np.sqrt(G_var+G_mean**2))
                sigma = np.sqrt(np.log(G_var/(G_mean**2)+1))
                self.G[l] = np.random.lognormal(mu, sigma,
                        size=(Nvol, Npart))

        # Define the model we're going to simulate
        self.m = modMPET("mpet", D=D)

    def SetUpParametersAndDomains(self):
        # Extract info from the config file
        # Simulation
        D = self.D
        # Geometry
        L_ref = D['L_ac'][1]
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
        csmax_ac = np.array([D['rhos_ac'][0], D['rhos_ac'][1]])/N_A
        # Ambipolar diffusivity
        Damb = ((zp+zm)*Dp*Dm)/(zp*Dp+zm*Dm)
        # Cation transference number
        tp = zp*Dp / (zp*Dp + zm*Dm)
        # Diffusive time scale
        td = L_ref**2 / Damb
        # Nondimensional Temperature
        T = float(D['Tabs'])/Tref
        # Electrode capacity ratio
        cap_c = (D['L_ac'][1] * (1-D['poros_ac'][1]) *
                D['P_L_ac'][1] * D['rhos_ac'][1])
        if D['Nvol_ac'][0]:
            # full porous anode with finite capacity
            cap_a = (D['L_ac'][0] * (1-D['poros_ac'][0]) *
                    D['P_L_ac'][0] * D['rhos_ac'][0])
            z = cap_c / cap_a
        else:
            # flat plate anode with assumed infinite supply of metal
            z = 0

        # Domains
        self.m.Ntrode.CreateArray(2)
        if self.Nvol_s >= 1:
            self.m.Nvol_s.CreateArray(self.Nvol_s)
        for l in self.trodes:
            self.m.Nvol_ac[l].CreateArray(self.Nvol_ac[l])
            self.m.Npart_ac[l].CreateArray(self.Npart_ac[l])
            for i in range(self.psd_num[l].shape[0]):
                for j in range(self.psd_num[l].shape[1]):
                    self.m.Nsld_mat_ac[l][i, j].CreateArray(
                            int(self.psd_num[l][i, j]))

        # Parameters
        self.m.T.SetValue(T)
        self.m.NumTrode.SetValue(len(self.trodes))
        for l in self.trodes:
            self.m.alpha_ac.SetValue(l, D['alpha_ac'][l])
            self.m.NumVol_ac.SetValue(l, self.Nvol_ac[l])
            self.m.NumPart_ac.SetValue(l, self.Npart_ac[l])
            self.m.L_ac.SetValue(l, D['L_ac'][l]/L_ref)
            self.m.dim_csmax_ac.SetValue(l, csmax_ac[l])
            self.m.poros_ac.SetValue(l, D['poros_ac'][l])
            self.m.epsbeta_ac.SetValue(l,
                    (1-D['poros_ac'][l]) * D['P_L_ac'][l] *
                    csmax_ac[l]/D['c0'])
            self.m.mcond_ac.SetValue(l,
                    D['mcond_ac'][l] * (td * k * N_A * Tref) /
                    (L_ref**2 * F**2 *D['c0']))
            if self.D['delPhiEqFit_ac'][l]:
                material = self.D['material_ac'][l]
                fits = delta_phi_fits.DPhiFits(self.D)
                phifunc = fits.materialData[material]
                self.m.dphi_eq_ref_ac.SetValue(l,
                        phifunc(self.D['cs0_ac'][l], 0))
            else:
                self.m.dphi_eq_ref_ac.SetValue(l,
                        0.0)
#                        -self.m.mu_reg_sln(D['cs0_ac'][l],
#                            D['Omga_ac'][l]/(k*Tref)))
            lmbda = D['lambda_ac'][l]/(k*Tref)
            self.m.lmbda_ac.SetValue(l, lmbda)
#            MHC_b = 2*np.sqrt(self.m.lmbda_ac.GetValue(l))
            MHC_erf_b = 2*np.sqrt(lmbda) # MHCerf
#            MHC_erf_b = np.sqrt(lmbda*np.pi) # MHCtanh
            self.m.MHC_erfstretch_ac.SetValue(l, MHC_erf_b)
#            self.m.MHC_Aa_ac.SetValue(l, )
            self.m.B_ac.SetValue(l,
                    D['B_ac'][l]/(k*Tref*D['rhos_ac'][l]))
            self.m.EvdW_ac.SetValue(l,
                    D['EvdW_ac'][l]/(k*Tref))
            self.m.cwet_ac.SetValue(l, D['cwet_ac'][l])
            for i in range(self.Nvol_ac[l]):
                for j in range(self.Npart_ac[l]):
                    solidShape = D['solidShape_ac'][l]
                    solidType = D['solidType_ac'][l]
                    p_num = float(self.psd_num[l][i, j])
                    p_len = self.psd_len[l][i, j]
                    # k0 is based on the _active_ area per volume for the region
                    # of the solid of interest.
                    if solidShape == "sphere":
                        # Spherical particles
                        p_area = (4*np.pi)*p_len**2
                        p_vol = (4./3)*np.pi*p_len**3
                    elif solidShape == "C3":
                        # C3 particles
                        p_area = 2 * 1.2263 * p_len**2
                        p_vol = 1.2263 * p_len**2 * D['partThick_ac'][l]
                    elif solidShape == "cylinder":
                        p_area = 2 * np.pi * p_len * D['partThick_ac'][l]
                        p_vol = np.pi * p_len**2 * D['partThick_ac'][l]
                    self.m.psd_num_ac[l].SetValue(i, j, p_num)
                    self.m.psd_len_ac[l].SetValue(i, j, p_len)
                    self.m.psd_area_ac[l].SetValue(i, j, p_area)
                    self.m.psd_vol_ac[l].SetValue(i, j, p_vol)
#                    kappatmp = D['kappa_ac'][l]/(k*Tref*D['rhos_ac'][l]*p_len**2)
#                    self.m.kappa_ac[l].SetValue(i, j, kappatmp)
                    self.m.kappa_ac[l].SetValue(i, j,
                            D['kappa_ac'][l]/(k*Tref*D['rhos_ac'][l]*p_len**2))
                    k0tmp = ((p_area/p_vol)*D['k0_ac'][l]*td)/(F*csmax_ac[l])
                    self.m.k0_ac[l].SetValue(i, j, k0tmp)
#                            ((p_area/p_vol)*D['k0_ac'][l]*td)/(F*csmax_ac[l]))
                    self.m.beta_s_ac[l].SetValue(i, j,
                            D['dgammasdc_ac'][l]*p_len*D['rhos_ac'][l]/D['kappa_ac'][l])
                    self.m.delta_L_ac[l].SetValue(i, j,
                            p_vol/(p_area*p_len))
                    self.m.MHC_Aa_ac[l].SetValue(i, j, # MHCerf
                            k0tmp / (spcl.erf(-lmbda/MHC_erf_b) + 1))
#                    self.m.MHC_Aa_ac[l].SetValue(i, j, # MHCtanh
#                            k0tmp / (np.tanh(-lmbda/MHC_erf_b) + 1))
                    self.m.scond_ac[l].SetValue(i, j,
                            D['scond_ac'][l] * (k*Tref)/(D['k0_ac'][l]*e*p_len**2))
                    self.m.Dsld_ac[l].SetValue(i, j,
                            D['Dsld_ac'][l]*td/p_len**2)
                    self.m.G_ac[l].SetValue(i, j,
                            self.G[l][i, j] * (k*Tref/e) * td/(F*csmax_ac[l]*p_vol))
                    if solidType in ["homog", "ACR", "CHR", "diffn"]:
                        self.m.Omga_ac[l].SetValue(i, j,
                                D['Omga_ac'][l]/(k*Tref))
                        self.m.Omgb_ac[l].SetValue(i, j,
                                D['Omgb_ac'][l]/(k*Tref))
                        self.m.Omgc_ac[l].SetValue(i, j,
                                D['Omgc_ac'][l]/(k*Tref))
                    elif solidType == "homog_sdn":
                        # Not sure about factor of nondimensional T. Thus,
                        # only use this when T = 1, Tabs = Tref = 298
                        self.m.Omga_ac[l].SetValue(i, j,
                                T*self.size2regsln(p_len))
                    else:
                        raise NotImplementedError("Solid types missing here")
        self.m.NumVol_s.SetValue(self.Nvol_s)
        self.m.L_s.SetValue(D['L_s']/L_ref)
        self.m.td.SetValue(td)
        self.m.zp.SetValue(zp)
        self.m.zm.SetValue(zm)
        self.m.tp.SetValue(tp)
        self.m.dim_Damb.SetValue(Damb)
        self.m.Dp.SetValue(Dp / Damb)
        self.m.Dm.SetValue(Dm / Damb)
        self.m.poros_s.SetValue(1.)
        self.m.phi_cathode.SetValue(0.)
        self.m.currset.SetValue(D['Crate']*td/3600)
        self.m.Vset.SetValue(D['Vset']*e/(k*Tref))
        self.m.z.SetValue(z)

    def SetUpVariables(self):
        Nlyte = self.Nvol_s + np.sum(self.Nvol_ac)
        phi_cathode = self.m.phi_cathode.GetValue()
        # Solids
        for l in self.trodes:
            cs0 = self.D['cs0_ac'][l]
            # Guess initial filling fractions
            self.m.ffrac_ac[l].SetInitialGuess(cs0)
            for i in range(self.Nvol_ac[l]):
                # Guess initial volumetric reaction rates
                self.m.j_plus_ac[l].SetInitialGuess(i, 0.0)
                # Guess initial value for the potential of the
                # electrodes
                if l == 1: # anode
                    self.m.phi_ac[l].SetInitialGuess(i, 0.0)
                else: # cathode
                    self.m.phi_ac[l].SetInitialGuess(i, phi_cathode)
                for j in range(self.Npart_ac[l]):
                    # Guess initial value for the average solid concentrations
                    self.m.c1bar_sld_ac[l].SetInitialGuess(i, j, cs0)
                    self.m.c2bar_sld_ac[l].SetInitialGuess(i, j, cs0)
                    # Set initial solid concentration values
                    Nij = self.m.Nsld_mat_ac[l][i, j].NumberOfPoints
                    epsrnd = 0.0001
                    rnd1 = epsrnd*(np.random.rand(Nij) - 0.5)
                    rnd2 = epsrnd*(np.random.rand(Nij) - 0.5)
                    rnd1 -= np.mean(rnd1)
                    rnd2 -= np.mean(rnd2)
                    for k in range(Nij):
                        self.m.c1_sld_ac[l][i, j].SetInitialCondition(k, cs0+rnd1[k])
                        self.m.c2_sld_ac[l][i, j].SetInitialCondition(k, cs0+rnd2[k])
        # Electrolyte
        c_lyte_init = self.D['c0']/1000. # normalize to 1 M = 1000 mol/m^3
        phi_guess = 0.
        for i in range(self.Nvol_s):
            self.m.c_lyte_s.SetInitialCondition(i, c_lyte_init)
            self.m.phi_lyte_s.SetInitialGuess(i, phi_guess)
        for l in self.trodes:
            for i in range(self.Nvol_ac[l]):
                self.m.c_lyte_ac[l].SetInitialCondition(i, c_lyte_init)
                self.m.phi_lyte_ac[l].SetInitialGuess(i, phi_guess)
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
        return param

    def test_input(self, D):
        if D['Tabs'] != 298 or D['Tref'] != 298:
            raise Exception("Temp dependence not implemented")
        if D['Nvol_ac'][1] < 1:
            raise Exception("Must have at least one porous electrode")
        for l in self.trodes:
            solidType = D['solidType_ac'][l]
            solidShape = D['solidShape_ac'][l]
            if D['simSurfCond_ac'][l] and solidType != "ACR":
                raise Exception("simSurfCond req. ACR")
            if solidType in ["ACR", "homog_sdn"] and solidShape != "C3":
                raise Exception("ACR and homog_sdn req. C3 shape")
            if solidType in ["CHR"] and solidShape not in ["sphere", "cylinder"]:
                raise NotImplementedError("CHR req. sphere or cylinder")
            if solidType not in ["ACR", "CHR", "homog", "homog_sdn", "diffn"]:
                raise NotImplementedError("Input solidType not defined")
            if solidShape not in ["C3", "sphere", "cylinder"]:
                raise NotImplementedError("Input solidShape not defined")
#            if solidType == "homog_sdn" and (D['Tabs'] != 298 or
#                    D['Tref'] != 298):
#                raise NotImplementedError("homog_snd req. Tref=Tabs=298")
            if solidType in ["diffn"] and solidShape != "sphere":
                raise NotImplementedError("diffn currently req. sphere")
            if D['delPhiEqFit_ac'][l] and solidType not in ["diffn", "homog"]:
                raise NotImplementedError("delPhiEqFit req. solidType = diffn or homog")
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
        self.interp = sint.interp1d(time_vec, noise_data, axis=0)
#        self.tlo = time_vec[0]
#        self.thi = time_vec[-1]
#        self.numnoise = len(time_vec)
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
        noise_vec = self.interp(time.Value)
        self.previous_output[:] = [time.Value, noise_vec] # it is a list now not a tuple
        self.counter += 1
        return adouble(noise_vec[self.position])
#        indx = (float(time.Value - self.tlo)/(self.thi-self.tlo) *
#                (self.numnoise - 1))
#        ilo = np.floor(indx)
#        ihi = np.ceil(indx)
#        # If we're exactly at a time in time_vec
#        if ilo == ihi:
#            noise_vec = self.noise_data[ilo, :]
#        else:
#            noise_vec = (self.noise_data[ilo, :] +
#                    (time.Value - self.time_vec[ilo]) /
#                    (self.time_vec[ihi] - self.time_vec[ilo]) *
#                    (self.noise_data[ihi, :] - self.noise_data[ilo, :])
#                    )
#        # previous_output is a reference to a common object and must
#        # be updated here - not deleted.  using self.previous_output = []
#        # it will delete the common object and create a new one
#        self.previous_output[:] = [time.Value, noise_vec] # it is a list now not a tuple
#        self.counter += 1
#        return adouble(float(noise_vec[self.position]))

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
            dkeybase = var.Name[len("mpet")+1:]
            mdict[dkeybase] = var.Values
            mdict[dkeybase + '_times'] = var.TimeValues
        try:
            import scipy.io
            scipy.io.savemat(self.ConnectionString,
                             mdict,
                             appendmat=False,
                             format='5',
                             long_field_names=False,
                             do_compression=False,
                             oned_as='row')
        except Exception, e:
            print 'Cannot call scipy.io.savemat(); is SciPy installed?\n' + str(e)

def setupDataReporters(simulation, outdir):
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
    # a hack to make compatible with pre/post r526 daetools
    try:
        simulation.dr.ConnectionString = simulation.dr.ConnectString
    except AttributeError:
        pass
    return datareporter

def consoleRun(D, outdir):
    # Create Log, Solver, DataReporter and Simulation object
    log          = daePythonStdOutLog()
    daesolver    = daeIDAS()
    simulation   = simMPET(D)
    datareporter = setupDataReporters(simulation, outdir)

    # Use SuperLU direct sparse LA solver
    lasolver = pySuperLU.daeCreateSuperLUSolver()
#    lasolver = pyTrilinos.daeCreateTrilinosSolver("Amesos_Umfpack", "")
    daesolver.SetLASolver(lasolver)
    
    # Enable reporting of all variables
    simulation.m.SetReportingOn(True)

    # Set relative tolerances
    daesolver.RelativeTolerance = 1e-6

    # Set the time horizon and the reporting interval
    # We need to get info about the system to figure out the
    # simulation time horizon
    Dp = D['Dp']
    Dm = D['Dm']
    zp = D['zp']
    zm = D['zm']
    Damb = ((zp+zm)*Dp*Dm)/(zp*Dp+zm*Dm)
    td = D['L_ac'][1]**2 / Damb
    currset = D['Crate'] * td/3600.
    if D['profileType'] == "CC" and currset != 0.0:
        simulation.TimeHorizon = np.abs(D['capFrac']/currset)
    else: # CV or zero current simulation
        simulation.TimeHorizon = D['tend']/td
    simulation.ReportingInterval = simulation.TimeHorizon/D['tsteps']
#    simulation.ReportingTimes = list(np.logspace(-4, np.log10(simulation.TimeHorizon), D['tsteps']))

    # Connect data reporter
    simName = simulation.m.Name + time.strftime(" [%d.%m.%Y %H:%M:%S]",
            time.localtime())
    if(datareporter.Connect("", simName) == False):
        sys.exit()

    # Initialize the simulation
    simulation.Initialize(daesolver, datareporter, log)

#    # Save model report
#    simulation.m.SaveModelReport(os.path.join(os.getcwd(),
#        "../2013_08_daetools/daetools_tutorials/report_mpet.xml"))

    # Solve at time=0 (initialization)
    simulation.SolveInitial()

    # Run
    try:
        simulation.Run()
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

    # Directories we'll store output in.
    outdir_name = "sim_output"
    outdir = os.path.join(os.getcwd(), outdir_name)
    # Make sure there's a place to store the output
    try:
        os.makedirs(outdir)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
        # Clean out the directory
        for file_obj in os.listdir(outdir):
            file_path = os.path.join(outdir, file_obj)
            if os.path.isfile(file_path):
                os.unlink(file_path)
            else:
                shutil.rmtree(file_path)
    paramFileName = "input_params.cfg"
    paramFile = os.path.join(outdir, paramFileName)
    IO.writeConfigFile(P, filename=paramFile)

    # Store info about this script
    try:
        # Git option, if it works -- commit info and current diff
        p1 = subprocess.Popen(['git', 'rev-parse', '--short', 'HEAD'],
                stdout=subprocess.PIPE)
        out1, err1 = p1.communicate()
        p2 = subprocess.Popen(['git', 'diff'],
                stdout=subprocess.PIPE)
        out2, err2 = p2.communicate()
        p3 = subprocess.Popen(['git', 'rev-parse', '--abbrev-ref', 'HEAD'],
                stdout=subprocess.PIPE)
        out3, err3 = p3.communicate()
        # Store commit info to file, as well as how to patch if
        # there's a diff
        with open(os.path.join(outdir, 'run_info.txt'), 'w') as fo:
            print >> fo, "branch name:"
            print >> fo, out3
            print >> fo, "commit hash:"
            print >> fo, out1
            print >> fo, "to run:"
            print >> fo, "$ git checkout [commit hash]"
            print >> fo, "$ patch -p1 < commit.diff:"
            print >> fo, "$ python[2] mpet.py input_params.cfg"
        with open(os.path.join(outdir, 'commit.diff'), 'w') as fo:
            print >> fo, out2
    except:
        # At least keep a copy of this file with the output
        shutil.copy(os.path.basename(__file__), outdir)
    try:
        shutil.copy("/etc/daetools/daetools.cfg", outdir)
    except:
        pass

    # Carry out the simulation
    consoleRun(D, outdir)

    # Final output for user
    if default_flag:
        print "\n\n*** WARNING: Used default file, ""{fname}"" ***".format(
                fname=default_file)
        print "Pass other parameter file as an argument to this script\n"
    else:
        print "\n\nUsed parameter file ""{fname}""\n\n".format(
                fname=paramfile)
    timeEnd = time.time()
    tTot = timeEnd - timeStart
    print "Total time:", tTot, "s"
    try:
        with open(os.path.join(outdir, 'run_info.txt'), 'a') as fo:
            print >> fo, "\nTotal run time:", tTot, "s"
    except Exception as e:
        pass

    # Copy simulation output to archive
    archivedir_name = time.strftime("%Y%m%d_%H%M%S", time.localtime())
    archivepath = os.path.join(os.getcwd(), "history")
    archivedir = os.path.join(archivepath, archivedir_name)
    try:
        os.makedirs(archivepath)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
    shutil.copytree(outdir, archivedir)

