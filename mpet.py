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

import mpet_params_IO
import delta_phi_fits
import mpetPorts
import mpetMaterials
import elyte_CST

eps = -1e-12

# Define some variable types
mole_frac_t = daeVariableType(name="mole_frac_t", units=unit(),
        lowerBound=0, upperBound=1, initialGuess=0.25,
        absTolerance=1e-6)
conc_t = daeVariableType(name="conc_t", units=unit(),
        lowerBound=0, upperBound=1e20, initialGuess=1.00,
        absTolerance=1e-6)
elec_pot_t = daeVariableType(name="elec_pot_t", units=unit(),
        lowerBound=-1e20, upperBound=1e20, initialGuess=0,
        absTolerance=1e-5)


class modMPET(daeModel):
    def __init__(self, Name, Parent=None, Description="", ndD_s=None,
            ndD_e=None):
        daeModel.__init__(self, Name, Parent, Description)

        if (ndD_s is None) or (ndD_e is None):
            raise Exception("Need input parameter dictionaries")
        self.ndD = ndD_s
        self.profileType = ndD_s['profileType']
        Nvol = ndD_s["Nvol"]
        Npart = ndD_s["Npart"]
        self.trodes = trodes = ndD_s["trodes"]

        # Domains where variables are distributed
        self.DmnCell = {} # domains over full cell dimensions
        self.DmnPart = {} # domains over particles in each cell volume
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

        # Variables
        self.c_lyte = {}
        self.phi_lyte = {}
        self.phi_bulk = {}
        self.phi_part = {}
        self.j_plus = {}
        self.ffrac = {}
        for l in trodes:
            # Concentration/potential in electrode regions of elyte
            self.c_lyte[l] = daeVariable("c_lyte_{l}".format(l=l),
                    conc_t, self,
                    "Concentration in the electrolyte in " +
                    "electrode {l}".format(l=l),
                    [self.DmnCell[l]])
            self.phi_lyte[l] = daeVariable("phi_lyte_{l}".format(l=l),
                    elec_pot_t, self,
                    "Electrostatic potential in electrolyte in " +
                    "electrode {l}".format(l=l),
                    [self.DmnCell[l]])
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
            self.c_lyte["s"] = daeVariable("c_lyte_s", conc_t, self,
                    "Concentration in the electrolyte in the separator",
                    [self.DmnCell["s"]])
            self.phi_lyte["s"] = daeVariable("phi_lyte_s", elec_pot_t, self,
                    "Electrostatic potential in electrolyte in separator",
                    [self.DmnCell["s"]])
        self.phi_applied = daeVariable("phi_applied", elec_pot_t, self,
                "Overall battery voltage (at anode current collector)")
        self.current = daeVariable("current", no_t, self,
                "Total current of the cell")
        self.dummyVar = daeVariable("dummyVar", no_t, self, "dummyVar")

        # Create models for representative particles within electrode
        # volumes and ports with which to talk to them.
        self.portsOutLyte = {}
        self.portsOutBulk = {}
        self.particles = {}
        for l in trodes:
            Nv = Nvol[l]
            Np = Npart[l]
            self.portsOutLyte[l] = np.empty(Nv, dtype=object)
            self.portsOutBulk[l] = np.empty((Nv, Np), dtype=object)
            self.particles[l] = np.empty((Nv, Np), dtype=object)
            for i in range(Nv):
                self.portsOutLyte[l][i] = mpetPorts.portFromElyte(
                        "portTrode{l}vol{i}".format(l=l,i=i),
                        eOutletPort, self,
                        "Electrolyte port to particles")
                for j in range(Np):
                    self.portsOutBulk[l][i, j] = mpetPorts.portFromBulk(
                        "portTrode{l}vol{i}part{j}".format(l=l,i=i,j=j),
                        eOutletPort, self,
                        "Bulk electrode port to particles")
                    solidType = ndD_e[l]["indvPart"][i, j]['type']
                    if solidType in ndD_s["2varTypes"]:
                        pMod = mpetMaterials.mod2var
                    elif solidType in ndD_s["1varTypes"]:
                        pMod = mpetMaterials.mod1var
                    else:
                        raise NotImplementedError("no model for given solid type")
                    self.particles[l][i, j] = pMod(
                            "partTrode{l}vol{i}part{j}".format(l=l,i=i,j=j),
                            self, ndD=ndD_e[l]["indvPart"][i, j],
                            ndD_s=ndD_s)
                    self.ConnectPorts(self.portsOutLyte[l][i],
                            self.particles[l][i, j].portInLyte)
                    self.ConnectPorts(
                            self.portsOutBulk[l][i, j],
                            self.particles[l][i, j].portInBulk)

    def DeclareEquations(self):
        daeModel.DeclareEquations(self)

        # Some values of domain lengths
        trodes = self.trodes
        ndD = self.ndD
        Nvol = ndD["Nvol"]
        Npart = ndD["Npart"]
        Nlyte = np.sum(Nvol.values())
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

        # Define the overall filling fraction in the electrodes
        for l in trodes:
            eq = self.CreateEquation("ffrac_{l}".format(l=l))
            eq.Residual = self.ffrac[l]()
            # Make a float of Vtot, total particle volume in electrode
            # Note: for some reason, even when "factored out", it's a bit
            # slower to use Sum(self.psd_vol_ac[l].array([], [])
            tmp = 0
            for i in range(Nvol[l]):
                for j in range(Npart[l]):
                    Vpart = ndD["psd_vol_FracTot"][l][i, j]
                    # For some reason the following slower, so
                    # it's helpful to factor out the Vtot sum:
                    # eq.Residual -= self.cbar_sld_ac[l](i, j) * (Vpart/Vtot)
                    tmp += self.particles[l][i, j].cbar() * Vpart
            eq.Residual -= tmp

        # Define dimensionless j_plus for each electrode volume
        for l in trodes:
            for i in range(Nvol[l]):
                eq = self.CreateEquation("j_plus_trode{l}vol{i}".format(i=i,l=l))
                # Start with no reaction, then add reactions for each
                # particle in the volume.
                res = 0
                # sum over particle volumes in given electrode volume
                for j in range(Npart[l]):
                    # The volume of this particular particle
                    Vj = ndD["psd_vol_FracVol"][l][i, j]
                    res += self.particles[l][i, j].dcbardt() * Vj
                eq.Residual = self.j_plus[l](i) - res

        # Define output port variables
        for l in trodes:
            for i in range(Nvol[l]):
                eq = self.CreateEquation("portout_c_trode{l}vol{i}".format(i=i,l=l))
                eq.Residual = (self.c_lyte[l](i) -
                        self.portsOutLyte[l][i].c_lyte())
                eq = self.CreateEquation("portout_p_trode{l}vol{i}".format(i=i,l=l))
                phi_lyte = self.phi_lyte[l](i)
                eq.Residual = (phi_lyte - self.portsOutLyte[l][i].phi_lyte())
                for j in range(Npart[l]):
                    eq = self.CreateEquation(
                            "portout_pm_trode{l}vol{i}part{j}".format(i=i,j=j,l=l))
                    eq.Residual = (self.phi_part[l](i, j) -
                            self.portsOutBulk[l][i, j].phi_m())

            # Simulate the potential drop along the bulk electrode
            # solid phase
            simBulkCond = ndD['simBulkCond'][l]
            if simBulkCond:
                # Calculate the RHS for electrode conductivity
                phi_tmp = np.empty(Nvol[l]+2, dtype=object)
                phi_tmp[1:-1] = [self.phi_bulk[l](i) for i in range(Nvol[l])]
                porosvec = np.empty(Nvol[l]+2, dtype=object)
                eps_sld = 1-self.ndD["poros"][l]
                porosvec[1:-1] = [eps_sld**(3./2) for i in range(Nvol[l])]
                porosvec[0] = porosvec[1]
                porosvec[-1] = porosvec[-2]
                porosvec = ((2*porosvec[1:]*porosvec[:-1])
                        / (porosvec[1:] + porosvec[:-1] + 1e-20))
                if l == "a": # anode
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
                RHS_phi_tmp = -np.diff(-porosvec*ndD["mcond"][l]*np.diff(phi_tmp)/dx)/dx
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
                    # charge conservation equation around this particle
                    eq = self.CreateEquation("phi_ac_trode{l}vol{i}part{j}".format(i=i,l=l,j=j))
                    if simPartCond:
                        # -dcsbar/dt = I_l - I_r
                        eq.Residual = (
                                self.particles[l][i, j].dcsbardt() + (
                                (-G_l * (phi_n - phi_l)) -
                                (-G_r * (phi_r - phi_n))))
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
                    if ndD["elyteModelType"] == "dilute":
                        eq.Residual = (ndD["poros"][l]*self.c_lyte[l].dt(i) +
                                ndD["epsbeta"][l]*(1-ndD["tp"])*self.j_plus[l](i) -
                                RHS_c[offset + i])
                    elif ndD["elyteModelType"] == "SM":
                        eq.Residual = (ndD["poros"][l]*self.c_lyte[l].dt(i) -
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
            eq.Residual = self.current() - ndD["currset"]*(
                    1 - np.exp(-Time()/(ndD["tend"]*1e-3)))
        elif self.profileType == "CV":
            # Keep applied potential constant
            eq = self.CreateEquation("applied_potential")
            eq.Residual = self.phi_applied() - ndD["Vset"]*(
                    1 - np.exp(-Time()/(ndD["tend"]*1e-3)))
#                    1)
#                    np.tanh(Time()/(45.0)))

        for eq in self.Equations:
            eq.CheckUnitsConsistency = False

        if self.profileType == "CC":
            # Set the condition to terminate the simulation upon reaching
            # a cutoff voltage.
            self.stopCondition = (
                    ((Abs(self.phi_applied()) <= ndD["phimin"]) |
                        (Abs(self.phi_applied()) >= ndD["phimax"]))
                    & (self.dummyVar() < 1))
            self.ON_CONDITION(self.stopCondition,
                    setVariableValues = [(self.dummyVar, 2)])

    def calc_lyte_RHS(self, cvec, phivec, Nvol, Nlyte):
        ndD = self.ndD
        zp = ndD["zp"]
        zm = ndD["zm"]
        nup = ndD["nup"]
        num = ndD["num"]
        nu = nup + num
        limtrode = ("c" if ndD["z"] < 1 else "a")
        # Discretization
        # The lengths are nondimensionalized by the cathode length
        dxvec = np.empty(np.sum(Nvol.values()) + 2, dtype=object)
        if Nvol["a"]:
            dxa = Nvol["a"] * [ndD["L"]["a"]/Nvol["a"]]
        else:
            dxa = []
        if Nvol["s"]:
            dxs = Nvol["s"] * [ndD["L"]["s"]/Nvol["s"]]
        else:
            dxs = []
        dxc = Nvol["c"] * [ndD["L"]["c"]/Nvol["c"]]
        dxtmp = np.array(dxa + dxs + dxc)
        dxvec[1:-1] = dxtmp
        dxvec[0] = dxvec[1]
        dxvec[-1] = dxvec[-2]
        dxd1 = (dxvec[0:-1] + dxvec[1:]) / 2.
        dxd2 = dxtmp

        # The porosity vector
        porosvec = np.empty(Nlyte + 2, dtype=object)
        # Use the Bruggeman relationship to approximate an effective
        # effect on the transport.
        porosvec[0:Nvol["a"]+1] = [ndD["poros"]["a"]**(3./2)
                for i in range(Nvol["a"]+1)] # anode
        porosvec[Nvol["a"]+1:Nvol["a"]+1 + Nvol["s"]] = [ndD["poros"]["s"]**(3./2)
                for i in range(Nvol["s"])] # separator
        porosvec[Nvol["a"]+1 + Nvol["s"]:] = [ndD["poros"]["c"]**(3./2)
                for i in range(Nvol["c"]+1)] # cathode
        poros_edges = (2*porosvec[1:]*porosvec[:-1])/(porosvec[1:] +
                porosvec[:-1] + 1e-20)

        if ndD["elyteModelType"] == "dilute":
            # Mass conservation equations
            ctmp = np.empty(Nlyte + 2, dtype=object)
            ctmp[1:-1] = cvec
            # If we don't have a real anode, the total current flowing
            # into the electrolyte is set
            if Nvol["a"] == 0:
                ctmp[0] = ( ctmp[1] + (self.current() *
                    ndD["epsbeta"][limtrode] *
                    (1-ndD["tp"])*dxvec[0])/poros_edges[0] )
            else: # porous anode -- no elyte flux at anode current collector
                ctmp[0] = ctmp[1]
            # No electrolyte flux at the cathode current collector
            ctmp[-1] = ctmp[-2]
            # Diffusive flux in the electrolyte
            cflux = -poros_edges*np.diff(ctmp)/dxd1
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
            c_edges = (2*ctmp[:-1]*ctmp[1:])/(ctmp[:-1] + ctmp[1:]+1e-20)
            Dp = ndD["Dp"]
            Dm = ndD["Dm"]
            # Typo in Todd's code in currdens equation
            currdens = (-((Dp - Dm)*np.diff(ctmp)/dxd1) -
                    (zp*Dp - zm*Dm)*c_edges*np.diff(phitmp)/dxd1)
            RHS_phi = -np.diff(poros_edges*currdens)/dxd2
            return (RHS_c, RHS_phi)

        elif ndD["elyteModelType"] == "SM":
            D, kappa, thermFac, tp0 = elyte_CST.getProps(ndD["SMset"])[:-1]
            # Vector of c values
            ctmp = np.empty(Nlyte + 2, dtype=object)
            ctmp[1:-1] = cvec
            # If we don't have a real anode, the total current flowing
            # into the electrolyte is set
            if Nvol["a"] == 0:
                ctmp[0] = ctmp[1] + (self.current() *
                    ndD["epsbeta"][limtrode] *
                    (1-tp0(ctmp[1]))) * (
                        dxvec[0]/(poros_edges[0]*D(ctmp[1])))
            else: # porous anode -- no elyte flux at anode current collector
                ctmp[0] = ctmp[1]
            # No electrolyte flux at the cathode current collector
            ctmp[-1] = ctmp[-2]

            # Vector of phi values and c at faces
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
            c_edges = (2*ctmp[:-1]*ctmp[1:])/(ctmp[:-1] + ctmp[1:]+1e-20)

            # current density in electrolyte
            i_edges = -poros_edges * kappa(c_edges) * (
                    np.diff(phitmp)/dxd1 -
                    nu/nup*(1-tp0(c_edges)) *
                    thermFac(c_edges) *
                    np.diff(np.log(ctmp))/dxd1
                    )
            # RHS for mass conservation equation
            RHS_c = (
                    np.diff(poros_edges * D(c_edges) *
                        np.diff(ctmp)/dxd1)/dxd2
                    +
                    1./(zp*nup)*np.diff((1-tp0(c_edges))*i_edges)/dxd2
                    )
            # RHS for charge conservation equation
            RHS_phi = -np.diff(i_edges)/dxd2
            return (RHS_c, RHS_phi)

class simMPET(daeSimulation):
    def __init__(self, ndD_s=None, ndD_e=None):
        daeSimulation.__init__(self)
        if (ndD_s is None) or (ndD_e is None):
            raise Exception("Need input parameter dictionaries")
        self.ndD_s = ndD_s
        self.ndD_e = ndD_e
        # Define the model we're going to simulate
        self.m = modMPET("mpet", ndD_s=ndD_s, ndD_e=ndD_e)

    def SetUpParametersAndDomains(self):
        # Domains
        ndD = self.ndD_s
        if ndD["Nvol"]["s"] >= 1:
            self.m.DmnCell["s"].CreateArray(ndD["Nvol"]["s"])
        for l in ndD["trodes"]:
            self.m.DmnCell[l].CreateArray(ndD["Nvol"][l])
            self.m.DmnPart[l].CreateArray(ndD["Npart"][l])
            for i in range(ndD["Nvol"][l]):
                for j in range(ndD["Npart"][l]):
                    self.m.particles[l][i, j].Dmn.CreateArray(
                            int(ndD["psd_num"][l][i, j]))

    def SetUpVariables(self):
        ndD_s = self.ndD_s
        Nvol = ndD_s["Nvol"]
        Npart = ndD_s["Npart"]
        Nlyte = np.sum(Nvol.values())
        phi_cathode = ndD_s["phi_cathode"]
        # Solids
        for l in ndD_s["trodes"]:
            cs0 = self.ndD_s['cs0'][l]
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
                    Nij = ndD_s["psd_num"][l][i, j]
                    # Guess initial value for the average solid concentrations
                    # and set initial value for solid concentrations
                    solidType = self.ndD_e[l]["indvPart"][i, j]["type"]
                    if solidType in ndD_s["1varTypes"]:
                        self.m.particles[l][i, j].cbar.SetInitialGuess(cs0)
                        for k in range(Nij):
                            self.m.particles[l][i, j].c.SetInitialCondition(k, cs0)
                    elif solidType in ndD_s["2varTypes"]:
                        self.m.particles[l][i, j].c1bar.SetInitialGuess(cs0)
                        self.m.particles[l][i, j].c2bar.SetInitialGuess(cs0)
                        self.m.particles[l][i, j].cbar.SetInitialGuess(cs0)
                        epsrnd = 0.0001
                        rnd1 = epsrnd*(np.random.rand(Nij) - 0.5)
                        rnd2 = epsrnd*(np.random.rand(Nij) - 0.5)
                        rnd1 -= np.mean(rnd1)
                        rnd2 -= np.mean(rnd2)
                        for k in range(Nij):
                            self.m.particles[l][i, j].c1.SetInitialCondition(k, cs0+rnd1[k])
                            self.m.particles[l][i, j].c2.SetInitialCondition(k, cs0+rnd2[k])
        # Electrolyte
        c_lyte_init = ndD_s['c0']
        phi_guess = 0.
        for i in range(Nvol["s"]):
            self.m.c_lyte["s"].SetInitialCondition(i, c_lyte_init)
            self.m.phi_lyte["s"].SetInitialGuess(i, phi_guess)
        for l in ndD_s["trodes"]:
            for i in range(Nvol[l]):
                self.m.c_lyte[l].SetInitialCondition(i, c_lyte_init)
                self.m.phi_lyte[l].SetInitialGuess(i, phi_guess)
        # Guess the initial cell voltage
        self.m.phi_applied.SetInitialGuess(0.0)
        self.m.dummyVar.AssignValue(0) # used for V cutoff condition

    def Run(self):
        """
        Overload the simulation "Run" function so that the simulation
        terminates when the specified condition is satisfied.
        """
        time = 0.
        while self.CurrentTime < self.TimeHorizon:
            nextTime = time + self.ReportingInterval
            if nextTime > self.TimeHorizon:
                nextTime = self.TimeHorizon
            self.Log.Message("Integrating from {t0:.2f} to {t1:.2f} s ...".format(
                t0=self.CurrentTime, t1=nextTime), 0)
            time = self.IntegrateUntilTime(nextTime,
                    eStopAtModelDiscontinuity, True)
            self.ReportData(self.CurrentTime)
            self.Log.SetProgress(int(100. * self.CurrentTime/self.TimeHorizon))
            if time < nextTime:
                # This means that the Integrate function returned
                # before reaching the specified nextTime.
                # This implies that the simulation stopped at the
                # "discontinuity".
                break

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

def consoleRun(ndD_s, ndD_e, outdir):
    # Create Log, Solver, DataReporter and Simulation object
    log          = daePythonStdOutLog()
    daesolver    = daeIDAS()
    simulation   = simMPET(ndD_s, ndD_e)
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
    simulation.TimeHorizon = ndD_s["tend"]
    simulation.ReportingInterval = ndD_s["tend"] / ndD_s['tsteps']

    # Connect data reporter
    simName = simulation.m.Name + time.strftime(" [%d.%m.%Y %H:%M:%S]",
            time.localtime())
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
    except KeyboardInterrupt:
        print "\nphi_applied at ctrl-C:", simulation.m.phi_applied.GetValue(), "\n"
        simulation.ReportData(simulation.CurrentTime)
    simulation.Finalize()

def main(paramfile="params_default.cfg", keepArchive=True):
    timeStart = time.time()
    # Get the parameters dictionary (and the config instance) from the
    # parameter file
    IO = mpet_params_IO.mpetIO()
    P_s, P_e = IO.getConfigs(paramfile)
    dD_s, ndD_s, dD_e, ndD_e = IO.getDictsFromConfigs(P_s, P_e)

    # Directories we'll store output in.
    outdir_name = time.strftime("%Y%m%d_%H%M%S", time.localtime())
    outdir_path = os.path.join(os.getcwd(), "history")
    outdir = os.path.join(outdir_path, outdir_name)
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
    paramFileName = "input_params_system.cfg"
    paramFile = os.path.join(outdir, paramFileName)
    IO.writeConfigFile(P_s, filename=paramFile)
    dictFile = os.path.join(outdir, "input_dict_system")
    IO.writeDicts(dD_s, ndD_s, filenamebase=dictFile)
    for trode in ndD_s["trodes"]:
        paramFileName = "input_params_{t}.cfg".format(t=trode)
        paramFile = os.path.join(outdir, paramFileName)
        IO.writeConfigFile(P_e[trode], filename=paramFile)
        dictFile = os.path.join(outdir,
                "input_dict_{t}".format(t=trode))
        IO.writeDicts(dD_e[trode], ndD_e[trode], filenamebase=dictFile)

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
    consoleRun(ndD_s, ndD_e, outdir)

    # Final output for user
    if paramfile == "params_default.cfg":
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

    # Copy simulation output to current directory
    tmpDir_name = "sim_output"
    tmpDir = os.path.join(os.getcwd(), tmpDir_name)
    try:
        shutil.rmtree(tmpDir)
    except OSError as exception:
        if exception.errno != errno.ENOENT:
            raise
    shutil.copytree(outdir, tmpDir)

    if not keepArchive:
        shutil.rmtree(outdir)

if __name__ == "__main__":
    default_file = "params_default.cfg"
    if len(sys.argv) < 2:
        paramfile = default_file
    else:
        paramfile = sys.argv[1]
    main(paramfile, keepArchive=True)
