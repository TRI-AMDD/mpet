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
        absTolerance=1e-6)
elec_pot_t = daeVariableType(name="elec_pot_t", units=unit(),
        lowerBound=-1e20, upperBound=1e20, initialGuess=0,
        absTolerance=1e-5)


class modMPET(daeModel):
    def __init__(self, Name, Parent=None, Description="", ndD=None):
        daeModel.__init__(self, Name, Parent, Description)

        if (ndD is None):
            raise Exception("Need input parameter dictionary")
        self.ndD = ndD
        self.profileType = ndD['profileType']
        Nvol = ndD["Nvol"]
        Npart = ndD["Npart"]
        self.trodes = trodes = ndD["trodes"]

        # Domains where variables are distributed
        self.DmnCell = {} # domains over full cell dimensions
        self.DmnPart = {} # domains over particles in each cell volume
        self.DmnPartSub = {} # domains within individual particles
        if Nvol["s"] >= 1: # If we have a separator
            self.DmnCell["s"] = daeDomain("DmnCell_s", self, unit(),
                    "Simulated volumes in the separator")
        for l in trodes:
            self.DmnCell[l] = daeDomain("DmnCell_{l}".format(l=l),
                    self, unit(),
                    "Simulated volumes in electrode " +
                    "{l}".format(l=l))
            self.DmnPart[l] = daeDomain("Npart_{l}".format(l=l),
                    self, unit(),
                    "Particles sampled in each control " +
                    "volume in electrode {l}".format(l=l))
            Nv = Nvol[l]
            Np = Npart[l]
            Nsld_mat = np.empty((Nv, Np), dtype=object)
            for i in range(Nv):
                for j in range(Np):
                    Nsld_mat[i, j] = daeDomain("trode{l}_vol{i}_part{j}".format(
                        i=i, j=j, l=l), self, unit(),
                        "Number of discretizations for particle "
                        + "j in volume i".format(i=i,j=j))
            self.DmnPartSub[l] = Nsld_mat

        # Variables
        self.c_lyte = {}
        self.phi_lyte = {}
        self.c_sld = {}
        self.phi_sld = {}
        self.cbar_sld = {}
        self.phi_bulk = {}
        self.phi_part = {}
        self.j_plus = {}
        self.ffrac = {}
        for l in trodes:
            # Concentration/potential in electrode regions of elyte
            self.c_lyte[l] = daeVariable("c_lyte_{l}".format(l=l),
                    mole_frac_t, self,
                    "Concentration in the electrolyte in " +
                    "electrode {l}".format(l=l),
                    [self.DmnCell[l]])
            self.phi_lyte[l] = daeVariable("phi_lyte_{l}".format(l=l),
                    elec_pot_t, self,
                    "Electrostatic potential in electrolyte in " +
                    "electrode {l}".format(l=l),
                    [self.DmnCell[l]])
            # Concentration in electrode active particles
            Nv = Nvol[l]
            Np = Npart[l]
            self.c_sld[l] = np.empty((Nv, Np), dtype=object)
            for i in range(Nv):
                for j in range(Np):
                    self.c_sld[l][i, j] = daeVariable(
                            "c_sld_trode{l}vol{i}part{j}".format(
                            i=i, j=j, l=l), mole_frac_t, self,
                            "Concentration in each solid particle",
                            [self.DmnPartSub[l][i, j]])
            # Potential in electrode active particles
            # Only make a variable of solid potentials if we have to
            # -- it's a lot of equations to keep track of for nothing
            # if we don't need it.
            if ndD['simSurfCond'][l]:
                self.phi_sld[l] = np.empty((Nv, Np), dtype=object)
                for i in range(Nv):
                    for j in range(Np):
                        self.phi_sld[l][i, j] = daeVariable(
                                "p_sld_trode{l}_vol{i}_part{j}".format(
                                i=i, j=j, l=l), elec_pot_t, self,
                                "Electrostatic potential in each solid particle",
                                [self.DmnPartSub[l][i, j]])
            else:
                self.phi_sld[l] = False
            # Average active particle concentrations
            self.cbar_sld[l] = daeVariable("cbar_sld_{l}".format(l=l),
                    mole_frac_t, self,
                    "Average concentration in each particle",
                    [self.DmnCell[l], self.DmnPart[l]])
            self.phi_bulk[l] = daeVariable("phi_bulk_{l}".format(l=l),
                    elec_pot_t, self,
                    "Electrostatic potential in the bulk solid",
                    [self.DmnCell[l]])
            self.phi_part[l] = daeVariable("phi_part_{l}".format(l=l),
                    elec_pot_t, self,
                    "Electrostatic potential at each particle",
                    [self.DmnCell[l], self.DmnPart[l]])
            self.j_plus[l] = daeVariable("j_plus_{l}".format(l=l),
                    no_t, self,
                    "Rate of reaction of positives per solid volume",
                    [self.DmnCell[l]])
            self.ffrac[l] = daeVariable("ffrac_{l}".format(l=l),
                mole_frac_t, self,
                "Overall filling fraction of solids in electrodes")
        if Nvol["s"] >= 1: # If we have a separator
            self.c_lyte["s"] = daeVariable("c_lyte_s", mole_frac_t, self,
                    "Concentration in the electrolyte in the separator",
                    [self.DmnCell["s"]])
            self.phi_lyte["s"] = daeVariable("phi_lyte_s", elec_pot_t, self,
                    "Electrostatic potential in electrolyte in separator",
                    [self.DmnCell["s"]])
        self.phi_applied = daeVariable("phi_applied", elec_pot_t, self,
                "Overall battery voltage (at anode current collector)")
        self.current = daeVariable("current", no_t, self,
                "Total current of the cell")

#        # Parameters
#        self.NumTrode = daeParameter("NumTrode", unit(), self,
#                "Number of electroes simulated (1 or 0)")
#        self.NumVol_ac = daeParameter("NumVol_ac", unit(), self,
#                "Number of volumes in the electrode",
#                [self.Ntrode])
#        self.NumPart_ac = daeParameter("NumPart_ac", unit(), self,
#                "Number of particles in each electrode volume",
#                [self.Ntrode])
#        self.NumVol_s = daeParameter("NumVol_s", unit(), self,
#                "Number of volumes in the electrolyte")
#        self.L_ac = daeParameter("L_ac", unit(), self,
#                "Length of electrodes (ndim to L_c)",
#                [self.Ntrode])
#        self.L_s = daeParameter("L_s", unit(), self,
#                "Length of separator (ndim to L_c)")
#        self.epsbeta_ac = daeParameter("epsbeta_ac", unit(), self,
#                "porosity times beta in electrodes",
#                [self.Ntrode])
#        self.zp = daeParameter("zp", unit(), self,
#                "cation charge number")
#        self.zm = daeParameter("zm", unit(), self,
#                "anion charge number")
#        self.tp = daeParameter("tp", unit(), self,
#                "positive transference number")
#        self.poros_s = daeParameter("poros_s", unit(), self,
#                "porosity in separator")
#        self.poros_ac = daeParameter("poros_ac", unit(), self,
#                "porosity in electrode",
#                [self.Ntrode])
#        self.phi_cathode = daeParameter("phi_cathode", unit(), self,
#                "reference potential, at the cathode + "
#                "(phi_applied is relative to this)")
#        self.td = daeParameter("td", unit(), self,
#                "Diffusive time [s]")
#        self.dim_Damb = daeParameter("dim_Damb", unit(), self,
#                "ambipolar diffusivity [m^2/s]")
#        self.Dp = daeParameter("Dp", unit(), self,
#                "non-dimensional diffusivity of positive ions")
#        self.Dm = daeParameter("Dm", unit(), self,
#                "non-dimensional diffusivity of negative ions")
#        self.Dsld_ac = np.empty(2, dtype=object)
#        self.kappa_ac = np.empty(2, dtype=object)
#        self.Omga_ac = np.empty(2, dtype=object)
#        self.k0_ac = np.empty(2, dtype=object)
#        self.beta_s_ac = np.empty(2, dtype=object)
#        self.delta_L_ac = np.empty(2, dtype=object)
#        self.MHC_Aa_ac = np.empty(2, dtype=object)
#        self.scond_ac = np.empty(2, dtype=object)
#        self.G_ac = np.empty(2, dtype=object)
#        self.psd_num_ac = np.empty(2, dtype=object)
#        self.psd_len_ac = np.empty(2, dtype=object)
#        self.psd_area_ac = np.empty(2, dtype=object)
#        self.psd_vol_ac = np.empty(2, dtype=object)
#        for l in trodes:
#            self.Dsld_ac[l] = daeParameter("Dsld_{l}".format(l=l),
#                    unit(), self,
#                    "Diffusivity in electrode active particles",
#                    [self.Nvol_ac[l], self.Npart_ac[l]])
#            self.kappa_ac[l] = daeParameter("kappa_{l}".format(l=l),
#                    unit(), self,
#                    "kappa for each particle",
#                    [self.Nvol_ac[l], self.Npart_ac[l]])
#            self.Omga_ac[l] = daeParameter("Omga_{l}".format(l=l),
#                    unit(), self,
#                    "regular solution parameter for each particle [J]",
#                    [self.Nvol_ac[l], self.Npart_ac[l]])
#            self.k0_ac[l] = daeParameter("k0_{l}".format(l=l),
#                    unit(), self,
#                    "exchange current density rate constant for each particle",
#                    [self.Nvol_ac[l], self.Npart_ac[l]])
#            self.beta_s_ac[l] = daeParameter("beta_s_{l}".format(l=l),
#                    unit(), self,
#                    "surface wetting, nondim: kappa*d(gamma_s)/dc",
#                    [self.Nvol_ac[l], self.Npart_ac[l]])
#            self.delta_L_ac[l] = daeParameter("delta_L_{l}".format(l=l),
#                    unit(), self,
#                    "Length ratios for particle: Vp/(Ap*Rp)",
#                    [self.Nvol_ac[l], self.Npart_ac[l]])
#            self.MHC_Aa_ac[l] = daeParameter("MHC_Aa_{l}".format(l=l),
#                    unit(), self,
#                    "MHC factor for erf approximation parameter",
##                    [self.Ntrode])
#                    [self.Nvol_ac[l], self.Npart_ac[l]])
#            self.scond_ac[l] = daeParameter("scond_{l}".format(l=l),
#                    unit(), self,
#                    "surface conductivity of particles",
#                    [self.Nvol_ac[l], self.Npart_ac[l]])
#            self.G_ac[l] = daeParameter("G_{l}".format(l=l),
#                    unit(), self,
#                    "conductance particles",
#                    [self.Nvol_ac[l], self.Npart_ac[l]])
#            self.psd_num_ac[l] = daeParameter("psd_numVols_{l}".format(l=l),
#                    unit(), self,
#                    "Particle numbers of discretizations",
#                    [self.Nvol_ac[l], self.Npart_ac[l]])
#            self.psd_len_ac[l] = daeParameter("psd_lengths_{l}".format(l=l),
#                    unit(), self,
#                    "Particle lengths [nm]",
#                    [self.Nvol_ac[l], self.Npart_ac[l]])
#            self.psd_area_ac[l] = daeParameter("psd_active_areas_{l}".format(l=l),
#                    unit(), self,
#                    "Particle active areas [nm^2]",
#                    [self.Nvol_ac[l], self.Npart_ac[l]])
#            self.psd_vol_ac[l] = daeParameter("psd_volumes_{l}".format(l=l),
#                    unit(), self,
#                    "Particle volumes [nm^3]",
#                    [self.Nvol_ac[l], self.Npart_ac[l]])
#        self.dphi_eq_ref_ac = daeParameter("dphi_eq_ref_ac",
#                unit(), self,
#                "dimensionless potential offset in referencing fit " +
#                "delta_phi_eq curves -- used only for initialization",
#                [self.Ntrode])
#        self.alpha_ac = daeParameter("alpha_ac", unit(), self,
#                " Charge transfer coefficient",
#                [self.Ntrode])
#        self.dim_csmax_ac = daeParameter("dim_csmax_ac", unit(), self,
#                "maximum lithium concentration in solid [mol/m^3]",
#                [self.Ntrode])
#        self.MHC_erfstretch_ac = daeParameter("MHC_b_ac", unit(), self,
#                "MHC factor for erf approximation parameter",
#                [self.Ntrode])
#        self.T = daeParameter("T", unit(), self,
#                "Non dimensional temperature")
#        self.currset = daeParameter("currset", unit(), self,
#                "dimensionless current")
#        self.Vset = daeParameter("Vset", unit(), self,
#                "dimensionless applied voltage (relative to " +
#                "Delta V OCV of the  cell)")
#        self.cwet_ac = daeParameter("cwet_ac", unit(), self,
#                "Wetted surface concentration",
#                [self.Ntrode])
#        self.B_ac = daeParameter("B_ac", unit(), self,
#                "Stress coefficient for each particle",
#                [self.Ntrode])
#        self.lmbda_ac = daeParameter("lmbda_ac", unit(), self,
#                "Marcus reorganizational energy",
#                [self.Ntrode])
#        self.mcond_ac = daeParameter("mcond_ac", unit(), self,
#                "conductivity of cathode",
#                [self.Ntrode])
#        self.z = daeParameter("z", unit(), self,
#                "capacity ratio of cathode:anode")

    def DeclareEquations(self):
        daeModel.DeclareEquations(self)

        # Some values of domain lengths
        trodes = self.trodes
#        Nvol_ac = [0, 0]
#        Npart_ac = [0, 0]
#        for l in trodes:
#            Nvol_ac[l] = self.Nvol_ac[l].NumberOfPoints
#            Npart_ac[l] = self.Npart_ac[l].NumberOfPoints
#        if self.D['Nvol_s'] >= 1: # if we have a separator
#            Nvol_s = self.Nvol_s.NumberOfPoints
#        else:
#            Nvol_s = 0
        ndD = self.ndD
        Nvol = ndD["Nvol"]
        Npart = ndD["Npart"]
        Nlyte = np.sum(Nvol.values())
        psd_num = ndD["psd_num"]
#        Nsld_mat_ac = np.empty(2, dtype=object)
#        for l in trodes:
#            Nsld_mat_ac[l] = np.zeros((Nvol_ac[l], Npart_ac[l]), dtype=np.integer)
#            for i in range(Nvol_ac[l]):
#                for j in range(Npart_ac[l]):
#                    Nptsij = self.Nsld_mat_ac[l][i, j].NumberOfPoints
#                    Nsld_mat_ac[l][i, j] = Nptsij
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
            for i in range(Nvol[l]):
                for j in range(Npart[l]):
                    Nij = psd_num[l][i, j]
                    eq = self.CreateEquation("cbar_trode{l}vol{i}part{j}".format(i=i,j=j,l=l))
                    r_vec, volfrac_vec = self.get_unit_solid_discr(
                            ndD['solidShape'][l],
                            ndD['solidType'][l], Nij)
                    eq.Residual = self.cbar_sld[l](i, j)
                    for k in range(Nij):
                        eq.Residual -= self.c_sld[l][i, j](k)*volfrac_vec[k]
#                    eq.BuildJacobianExpressions = True

        # Define the overall filling fraction in the electrodes
        for l in trodes:
            eq = self.CreateEquation("ffrac_{l}".format(l=l))
            eq.Residual = self.ffrac[l]()
            # Make a float of Vtot, total particle volume in electrode
            # Note: for some reason, even when "factored out", it's a bit
            # slower to use Sum(self.psd_vol_ac[l].array([], [])
#            Vtot = np.sum(psd_vol[l])
            tmp = 0
            for i in range(Nvol[l]):
                for j in range(Npart[l]):
#                    Vpart = psd_vol[l][i, j]
                    Vpart = ndD["psd_vol_FracTot"][l][i, j]
                    # For some reason the following slower, so
                    # it's helpful to factor out the Vtot sum:
                    # eq.Residual -= self.cbar_sld_ac[l](i, j) * (Vpart/Vtot)
                    tmp += self.cbar_sld[l](i, j) * Vpart
#            eq.Residual -= tmp / Vtot
            eq.Residual -= tmp

        # Define dimensionless j_plus for each electrode volume
        for l in trodes:
            for i in range(Nvol[l]):
                eq = self.CreateEquation("j_plus_trode{l}vol{i}".format(i=i,l=l))
                # Start with no reaction, then add reactions for each
                # particle in the volume.
                res = 0
                # sum over particle volumes in given electrode volume
#                Vu = Sum(self.psd_vol_ac[l].array(i, []))
#                Vu = np.sum(psd_vol[l][i, :])
                for j in range(Npart[l]):
                    # The volume of this particular particle
#                    Vj = psd_vol[l][i, j]
                    Vj = ndD["psd_vol_FracVol"][l][i, j]
                    Nij = psd_num[l][i, j]
                    r_vec, volfrac_vec = self.get_unit_solid_discr(
                            ndD['solidShape'][l],
                            ndD['solidType'][l], Nij)
                    tmp = 0
                    for k in range(Nij):
                        tmp += (self.c_sld[l][i, j].dt(k) *
                                volfrac_vec[k])
                    res += tmp * Vj
#                eq.Residual = self.j_plus[l](i) - res/Vu
                eq.Residual = self.j_plus[l](i) - res

        # Solid active particle concentrations, potential, and bulk
        # solid potential
        for l in trodes:
            # Solid active particles concentration/potential
            for i in range(Nvol[l]):
                # Calculate the solid concentration rates of change
                for j in range(Npart[l]):
                    # Prepare the RHS function
                    Nij = psd_num[l][i, j]
                    # Note Mmat is often the identity matrix...
                    (Mmat, RHS_c_sld_ij) = self.calc_sld_dcs_dt(l, i, j)
                    dcdt_vec = np.empty(Nij, dtype=object)
                    dcdt_vec[0:Nij] = [self.c_sld[l][i, j].dt(k) for k in range(Nij)]
                    LHS_vec = self.MX(Mmat, dcdt_vec)
#                    noisevec = self.noise_local[l][i, j] # NOISE
                    # Set up equations: Mmat*dcdt = RHS
                    for k in range(Nij):
                        eq = self.CreateEquation(
                                "dcsdt_trode{l}vol{i}part{j}discr{k}".format(
                                    i=i,j=j,k=k,l=l))
#                        eq.Residual = self.c_sld[i, j].dt(k) - RHS_c_sld_ij[k]
                        eq.Residual = LHS_vec[k] - RHS_c_sld_ij[k]
#                        eq.Residual = (LHS_vec[k] - RHS_c_sld_ij[k] +
#                                noisevec[k]()) # NOISE

                # Also calculate the potential drop along cathode
                # particle surfaces, if desired
                simSurfCathCond = ndD['simSurfCond'][l]
                if simSurfCathCond:
                    # Conservation of charge in the solid particles with
                    # Ohm's Law
                    LHS = self.calc_part_surf_LHS(l, i, j)
                    k0_part = ndD["k0"][l][i, j]
                    for k in range(Nij):
                        eq = self.CreateEquation(
                                "charge_cons_trode{l}vol{i}part{j}discr{k}".format(
                                    i=i,j=j,k=k,l=l))
                        RHS = self.c_sld[l][i, j].dt(k) / k0_part
                        eq.Residual = LHS[k] - RHS

            # Simulate the potential drop along the bulk electrode
            # solid phase
            simBulkCond = ndD['simBulkCond'][l]
            if simBulkCond:
                # Calculate the RHS for electrode conductivity
                phi_tmp = np.empty(Nvol[l]+2, dtype=object)
                phi_tmp[1:-1] = [self.phi_bulk[l](i) for i in range(Nvol[l])]
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
                    phi_tmp[-1] = ndD["phi_cathode"]
                dx = 1./Nvol[l]
                RHS_phi_tmp = -np.diff(-ndD["mcond"][l]*np.diff(phi_tmp)/dx)/dx
            # Actually set up the equations for bulk solid phi
            for i in range(Nvol[l]):
                eq = self.CreateEquation("phi_ac_trode{l}vol{i}".format(i=i,l=l))
                if simBulkCond:
                    eq.Residual = (-ndD["epsbeta"][l]*self.j_plus[l](i) -
                            RHS_phi_tmp[i])
                else:
                    if l == "a": # anode
                        eq.Residual = self.phi_bulk[l](i) - self.phi_applied()
                    else: # cathode
                        eq.Residual = self.phi_bulk[l](i) - ndD["phi_cathode"]

            # Simulate the potential drop along the connected
            # particles
            simPartCond = ndD['simPartCond'][l]
            for i in range(Nvol[l]):
                phi_bulk = self.phi_bulk[l](i)
                for j in range(Npart[l]):
                    G_l = ndD["G"][l][i, j]
                    phi_n = self.phi_part[l](i, j)
                    if j == 0: # reference bulk phi
                        phi_l = self.phi_bulk[l](i)
                    else:
                        phi_l = self.phi_part[l](i, j - 1)
                    if j == (Npart[l] - 1): # No particle at and of "chain"
                        G_r = 0
                        phi_r = phi_n
                    else:
                        G_r = ndD["G"][l][i, j + 1]
                        phi_r = self.phi_part[l](i, j + 1)
                    # Find average dcs/dt for this particle
                    Nij = psd_num[l][i, j]
                    r_vec, volfrac_vec = self.get_unit_solid_discr(
                            ndD['solidShape'][l],
                            ndD['solidType'][l], Nij)
                    dcsbardt = 0
                    for k in range(Nij):
                        dcsbardt += (self.c_sld[l][i, j].dt(k) *
                                volfrac_vec[k])
                    # charge conservation equation around this particle
                    eq = self.CreateEquation("phi_ac_trode{l}vol{i}part{j}".format(i=i,l=l,j=j))
                    if simPartCond:
                        # -dcsbar/dt = I_l - I_r
                        eq.Residual = dcsbardt + (
                                (-G_l * (phi_n - phi_l)) -
                                (-G_r * (phi_r - phi_n)))
                    else:
                        eq.Residual = self.phi_part[l](i, j) - self.phi_bulk[l](i)

        # If we have a single electrode volume (in a perfect bath),
        # electrolyte equations are simple
        if Nvol["a"] == 0 and Nvol["s"] == 0 and Nvol["c"] == 1:
            eq = self.CreateEquation("c_lyte")
            eq.Residual = self.c_lyte["c"].dt(0) - 0
            eq = self.CreateEquation("phi_lyte")
            eq.Residual = self.phi_lyte["c"](0) - self.phi_applied()
        else:
            # Calculate RHS for electrolyte equations
            c_lyte = np.empty(Nlyte, dtype=object)
            c_lyte[0:Nvol["a"]] = [self.c_lyte["a"](i)
                    for i in range(Nvol["a"])] # anode
            c_lyte[Nvol["a"]:Nvol["a"] + Nvol["s"]] = [self.c_lyte["s"](i)
                    for i in range(Nvol["s"])] # separator
            c_lyte[Nvol["a"] + Nvol["s"]:Nlyte] = [self.c_lyte["c"](i)
                    for i in range(Nvol["c"])] # cathode
            phi_lyte = np.empty(Nlyte, dtype=object)
            phi_lyte[0:Nvol["a"]] = [self.phi_lyte["a"](i)
                    for i in range(Nvol["a"])] # anode
            phi_lyte[Nvol["a"]:Nvol["a"] + Nvol["s"]] = [self.phi_lyte["s"](i)
                    for i in range(Nvol["s"])] # separator
            phi_lyte[Nvol["a"] + Nvol["s"]:Nlyte] = [self.phi_lyte["c"](i)
                    for i in range(Nvol["c"])] # cathode
            (RHS_c, RHS_phi) = self.calc_lyte_RHS(c_lyte, phi_lyte,
                    Nvol, Nlyte)
            # Equations governing the electrolyte in the separator
            offset = Nvol["a"]
            for i in range(Nvol["s"]):
                # Mass Conservation
                eq = self.CreateEquation(
                        "sep_lyte_mass_cons_vol{i}".format(i=i))
                eq.Residual = (ndD["poros"]["s"]*self.c_lyte["s"].dt(i) -
                        RHS_c[offset + i])
                # Charge Conservation
                eq = self.CreateEquation(
                        "sep_lyte_charge_cons_vol{i}".format(i=i))
                eq.Residual = (RHS_phi[offset + i])
            # Equations governing the electrolyte in the electrodes.
            # Here, we are coupled to the total reaction rates in the
            # solids.
            for l in trodes:
                if l == "a": # anode
                    offset = 0
                else: # cathode
                    offset = Nvol["a"] + Nvol["s"]
                for i in range(Nvol[l]):
                    # Mass Conservation
                    eq = self.CreateEquation(
                            "lyteMassCons_trode{l}vol{i}".format(i=i,l=l))
                    eq.Residual = (ndD["poros"][l]*self.c_lyte[l].dt(i) +
                            ndD["epsbeta"][l]*(1-ndD["tp"])*self.j_plus[l](i) -
                            RHS_c[offset + i])
                    # Charge Conservation
                    eq = self.CreateEquation(
                            "lyteChargeCons_trode{l}vol{i}".format(i=i,l=l))
                    eq.Residual = (ndD["epsbeta"][l]*self.j_plus[l](i) -
                            RHS_phi[offset + i])

        # Define the total current. This can be done in either anode
        # or cathode equivalently.
        eq = self.CreateEquation("Total_Current")
        eq.Residual = self.current()
        limtrode = ("c" if ndD["z"] < 1 else "a")
        dx = 1./Nvol[limtrode]
        for i in range(Nvol[limtrode]):
            if limtrode == "a":
                eq.Residual += dx * self.j_plus[limtrode](i)
            else:
                eq.Residual -= dx * self.j_plus[limtrode](i)

        if self.profileType == "CC":
            # Total Current Constraint Equation
            eq = self.CreateEquation("Total_Current_Constraint")
#            if ndD["currset"] != 0.0:
#                timeHorizon = 1./np.abs(ndD["currset"])
#            else:
#                timeHorizon = ndD["tend"]
            eq.Residual = self.current() - ndD["currset"]*(
                    1 - np.exp(-Time()/(ndD["tend"]*1e-3)))
        elif self.profileType == "CV":
            # Keep applied potential constant
            eq = self.CreateEquation("applied_potential")
#            timeHorizon = ndD["tend"]
            eq.Residual = self.phi_applied() - ndD["Vset"]*(
                    1 - np.exp(-Time()/(ndD["tend"]*1e-3)))

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
        ndD = self.ndD
        simSurfCond = ndD['simSurfCond'][l]
        solidType = ndD['solidType'][l]
        solidShape = ndD['solidShape'][l]
        rxnType = ndD['rxnType'][l]
        delPhiEqFit = ndD['delPhiEqFit'][l]
        # Get variables for this particle/electrode volume
        phi_lyte = self.phi_lyte[l](i)
        phi_m = self.phi_part[l](i, j)
        c_lyte = self.c_lyte[l](i)
        # Get the relevant parameters for this particle
        cbar = self.cbar_sld[l](i, j) # only used for ACR/CHR
        k0 = ndD["k0"][l][i, j]
        kappa = ndD["kappa"][l][i, j] # only used for ACR/CHR
        lmbda = ndD["lambda"][l] # Only used for Marcus/MHC
        Aa = ndD["MHC_Aa"][l][i, j] # Only used for MHC
        b = ndD["MHC_erfstretch"][l] # Only used for MHC
        alpha = ndD["alpha"][l] # Only used for BV
        Omga = ndD["Omga"][l][i, j]
        Ds = ndD["Dsld"][l][i, j] # Only used for "diffn"
        # We need the (non-dimensional) temperature to get the
        # reaction rate dependence correct
        T = ndD["T"]
        # Number of volumes in current particle
        Nij = ndD["psd_num"][l][i, j]
        # Concentration (profile?) in the solid
        c_sld = np.empty(Nij, dtype=object)
        c_sld[:] = [self.c_sld[l][i, j](k) for k in range(Nij)]
        # Assume dilute electrolyte
        act_O = c_lyte
        mu_O = T*np.log(Max(eps, act_O))
        # Calculate chemical potential of reduced state

        if solidType in ["ACR", "homog", "homog_sdn"]:
            if solidType == "ACR":
                # Make a blank array to allow for boundary conditions
                cstmp = np.empty(Nij+2, dtype=object)
                cstmp[1:-1] = c_sld
                cstmp[0] = ndD["cwet"][l]
                cstmp[-1] = ndD["cwet"][l]
                dxs = 1./Nij
                curv = np.diff(cstmp, 2)/(dxs**2)
                mu_R = ( self.mu_reg_sln(c_sld, Omga) - kappa*curv
                        + ndD["B"][l]*(c_sld - cbar) )
                # If we're also simulating potential drop along the solid,
                # use that instead of self.phi_c(i)
                if simSurfCond:
                    phi_m = np.empty(Nij, dtype=object)
                    phi_m[:] = [self.phi_sld[l][i, j](k) for k in range(Nij)]
            elif solidType in ["homog", "homog_sdn"]:
                mu_R = self.mu_reg_sln(c_sld, Omga)
            # XXX -- Temp dependence!
            act_R = np.exp(mu_R/T)
            # eta = electrochem pot_R - electrochem pot_O
            # eta = (mu_R + phi_R) - (mu_O + phi_O)
            if delPhiEqFit:
                material = ndD['material'][l]
                fits = delta_phi_fits.DPhiFits(ndD["T"])
                phifunc = fits.materialData[material]
                delta_phi_eq = phifunc(c_sld[-1], ndD["dphi_eq_ref"][l])
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
                RHS = np.empty(Nij, dtype=object)
                c_diffs = np.diff(c_sld)
                RHS[1:Nij - 1] = 4*np.pi*(
                        Ds*(r_vec[1:Nij - 1] + dr/2.)**2*c_diffs[1:]/dr -
                        Ds*(r_vec[1:Nij - 1] - dr/2.)**2*c_diffs[:-1]/dr )
                RHS[0] = 4*np.pi*Ds*(dr/2.)**2*c_diffs[0]/dr
                # Take the surface concentration
                c_surf = c_sld[-1]
                # Assuming dilute solid
                act_R_surf = c_surf
                mu_R_surf = T*np.log(Max(eps, act_R_surf))
            # CHR
            elif solidType in ["CHR"]:
                # mu_R is like for ACR, except kappa*curv term
                mu_R = ( self.mu_reg_sln(c_sld, Omga) +
                        ndD["B"][l]*(c_sld - cbar) )
                # Surface conc gradient given by natural BC
                beta_s = ndD["beta_s"][l][i, j]
                if solidShape == "sphere":
                    mu_R[0] -= 3 * kappa * (2*c_sld[1] - 2*c_sld[0])/dr**2
                    mu_R[1:Nij - 1] -= kappa * (np.diff(c_sld, 2)/dr**2 +
                            (c_sld[2:] - c_sld[0:-2])/(dr*r_vec[1:-1])
                            )
                    mu_R[Nij - 1] -= kappa * ((2./Rs)*beta_s +
                            (2*c_sld[-2] - 2*c_sld[-1] + 2*dr*beta_s)/dr**2
                            )
                elif solidShape == "cylinder":
                    mu_R[0] -= 2 * kappa * (2*c_sld[1] - 2*c_sld[0])/dr**2
                    mu_R[1:Nij - 1] -= kappa * (np.diff(c_sld, 2)/dr**2 +
                            (c_sld[2:] - c_sld[0:-2])/(2 * dr*r_vec[1:-1])
                            )
                    mu_R[Nij - 1] -= kappa * ((1./Rs)*beta_s +
                            (2*c_sld[-2] - 2*c_sld[-1] + 2*dr*beta_s)/dr**2
                            )
                mu_R_surf = mu_R[-1]
                # With chem potentials, can now calculate fluxes
                Flux_vec = np.empty(Nij+1, dtype=object)
                Flux_vec[0] = 0  # Symmetry at r=0
                c_edges = (c_sld[0:-1]  + c_sld[1:])/2.
                Flux_vec[1:Nij] = (Ds/T * (1 - c_edges) * c_edges *
                        np.diff(mu_R)/dr)
                # Take the surface concentration
                c_surf = c_sld[-1]
                act_R_surf = np.exp(mu_R_surf/T)
            # Figure out overpotential
            if delPhiEqFit:
                material = ndD['material'][l]
                fits = delta_phi_fits.DPhiFits(ndD["T"])
                phifunc = fits.materialData[material]
                delta_phi_eq = phifunc(c_surf, ndD["dphi_eq_ref"][l])
                eta = (phi_m - phi_lyte) - delta_phi_eq
            else:
                eta = (mu_R_surf + phi_m) - (mu_O + phi_lyte)
            # Calculate reaction rate
            if rxnType == "Marcus":
                Rxn = self.R_Marcus(k0, lmbda, c_lyte, c_surf, eta, T)
            elif rxnType == "BV":
                Rxn = self.R_BV(k0, alpha, c_surf, act_O, act_R_surf, eta, T)
            elif rxnType == "MHC":
                Rxn = self.R_MHC(k0, lmbda, eta, Aa, b, T, c_surf)
            # Finish up RHS discretization at particle surface
            if solidType in ["diffn"]:
                RHS[-1] = 4*np.pi*(Rs**2 * ndD["delta_L"][l][i, j] * Rxn -
                        Ds*(Rs - dr/2)**2*c_diffs[-1]/dr )
            elif solidType in ["CHR"]:
                Flux_vec[Nij] = ndD["delta_L"][l][i, j] * Rxn
                if solidShape == "sphere":
                    area_vec = 4*np.pi*edges**2
                elif solidShape == "cylinder":
                    area_vec = 2*np.pi*edges  # per unit height
                RHS = np.diff(Flux_vec*area_vec)
            return (M, RHS)

    def calc_lyte_RHS(self, cvec, phivec, Nvol, Nlyte):
        # Discretization
        # The lengths are nondimensionalized by the cathode length
        dxvec = np.empty(np.sum(Nvol.values()) + 2, dtype=object)
        dxa = np.array(Nvol["a"] * [self.ndD["L"]["a"]/Nvol["a"]])
        dxc = np.array(Nvol["c"] * [self.ndD["L"]["c"]/Nvol["c"]])
        dxs = np.array(Nvol["s"] * [self.ndD["L"]["s"]/Nvol["s"]])
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
        porosvec[0:Nvol["a"]] = [self.ndD["poros"]["a"]**(3./2)
                for i in range(Nvol["a"])] # anode
        porosvec[Nvol["a"]:Nvol["a"] + Nvol["s"]] = [self.ndD["poros"]["s"]**(3./2)
                for i in range(Nvol["s"])] # separator
        porosvec[Nvol["a"] + Nvol["s"]:Nlyte+1] = [self.ndD["poros"]["c"]**(3./2)
                for i in range(Nvol["c"]+1)] # cathode

        # Mass conservation equations
        ctmp = np.empty(Nlyte + 2, dtype=object)
        ctmp[1:-1] = cvec
        # If we don't have a real anode, the total current flowing
        # into the electrolyte is set
        limtrode = ("c" if self.ndD["z"] < 1 else "a")
        if Nvol["a"] == 0:
            ctmp[0] = ( ctmp[1] + (self.current() *
                self.ndD["epsbeta"][limtrode] *
                (1-self.ndD["tp"])*dxvec[0])/porosvec[0] )
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
        if Nvol["a"] == 0:
            phitmp[0] = self.phi_applied()
        else: # porous anode -- no flux into anode current collector
            phitmp[0] = phitmp[1]
        # No flux into cathode current collector from the electrolyte
        phitmp[-1] = phitmp[-2]
        # We need average values of c_lyte for the current densities
        # at the finite volume boundaries
        c_edges = (ctmp[0:-1] + ctmp[1:])/2.
        zp = self.ndD["zp"]
        zm = self.ndD["zm"]
        Dp = self.ndD["Dp"]
        Dm = self.ndD["Dm"]
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
        Nij = ndD["psd_num"][l][i, j]
        # solid potential variables for this particle
        phi_tmp = np.empty(Nij + 2, dtype=object)
        phi_tmp[1:-1] = [self.phi_sld[l][i, j](k) for k in
                range(Nij)]
        # BC's -- "touching carbon at each end"
        phi_s_local = self.phi_part[l](i, j)
        phi_tmp[0] = phi_s_local
        phi_tmp[-1] = phi_s_local
        # LHS
        dx = 1./Nij
        phi_edges = (phi_tmp[0:-1] + phi_tmp[1:])/2.
#        curr_dens = -self.scond(i, j)*np.diff(phi_tmp, 1)/dx
        # XXX -- Temp dependence!
        scond_vec = self.ndD["scond"][l][i, j]*np.exp(-1*(phi_edges -
                phi_s_local))
        curr_dens = -scond_vec*np.diff(phi_tmp, 1)/dx
        return np.diff(curr_dens, 1)/dx

    def mu_reg_sln(self, c, Omga):
        return np.array([ Omga*(1-2*c[i])
                + self.ndD["T"]*Log(Max(eps, c[i])/Max(eps, 1-c[i]))
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
    def __init__(self, ndD=None):
        daeSimulation.__init__(self)
        if ndD is None:
            raise Exception("Need parameters input")
        self.ndD = ndD
        # Define the model we're going to simulate
        self.m = modMPET("mpet", ndD=ndD)

    def SetUpParametersAndDomains(self):
        # Domains
        ndD = self.ndD
        if ndD["Nvol"]["s"] >= 1:
            self.m.DmnCell["s"].CreateArray(ndD["Nvol"]["s"])
        for l in ndD["trodes"]:
            self.m.DmnCell[l].CreateArray(ndD["Nvol"][l])
            self.m.DmnPart[l].CreateArray(ndD["Npart"][l])
            for i in range(ndD["Nvol"][l]):
                for j in range(ndD["Npart"][l]):
                    self.m.DmnPartSub[l][i, j].CreateArray(
                            int(ndD["psd_num"][l][i, j]))

    def SetUpVariables(self):
        ndD = self.ndD
        Nvol = ndD["Nvol"]
        Npart = ndD["Npart"]
        Nlyte = np.sum(Nvol.values())
        phi_cathode = ndD["phi_cathode"]
        # Solids
        for l in ndD["trodes"]:
            cs0 = ndD['cs0'][l]
            # Guess initial filling fractions
            self.m.ffrac[l].SetInitialGuess(cs0)
            for i in range(Nvol[l]):
                # Guess initial volumetric reaction rates
                self.m.j_plus[l].SetInitialGuess(i, 0.0)
                # Guess initial value for the potential of the
                # electrodes
                if l == "a": # anode
                    self.m.phi_bulk[l].SetInitialGuess(i, 0.0)
                else: # cathode
                    self.m.phi_bulk[l].SetInitialGuess(i, phi_cathode)
                for j in range(Npart[l]):
                    # Guess initial value for the average solid concentrations
                    self.m.cbar_sld[l].SetInitialGuess(i, j, cs0)
                    # Set initial solid concentration values
                    Nij = ndD["psd_num"][l][i, j]
                    for k in range(Nij):
                        self.m.c_sld[l][i, j].SetInitialCondition(k, cs0)
        # Electrolyte
        c_lyte_init = ndD['c0']
        phi_guess = 0.
        for i in range(Nvol["s"]):
            self.m.c_lyte["s"].SetInitialCondition(i, c_lyte_init)
            self.m.phi_lyte["s"].SetInitialGuess(i, phi_guess)
        for l in ndD["trodes"]:
            for i in range(Nvol[l]):
                self.m.c_lyte[l].SetInitialCondition(i, c_lyte_init)
                self.m.phi_lyte[l].SetInitialGuess(i, phi_guess)
        # Guess the initial cell voltage
        self.m.phi_applied.SetInitialGuess(0.0)

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

def consoleRun(ndD, outdir):
    # Create Log, Solver, DataReporter and Simulation object
    log          = daePythonStdOutLog()
    daesolver    = daeIDAS()
    simulation   = simMPET(ndD)
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
#    if ndD['profileType'] == "CC" and ndD["currset"] != 0.0:
#        simulation.TimeHorizon = 1./np.abs(ndD['currset'])
#    else: # CV or zero current simulation
#        simulation.TimeHorizon = ndD['tend']
    simulation.TimeHorizon = ndD["tend"]
    simulation.ReportingInterval = simulation.TimeHorizon/ndD['tsteps']

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
    dD, ndD = IO.getDictFromConfig(P)

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
    dictFile = os.path.join(outdir, "input_dict")
    IO.writeDicts(dD, ndD, filenamebase=dictFile)

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
    consoleRun(ndD, outdir)

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

