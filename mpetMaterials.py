#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import numpy as np
import scipy.sparse as sprs
import scipy.special as spcl

from daetools.pyDAE import *
from daetools.pyDAE.data_reporters import *
from daetools.solvers.superlu import pySuperLU

import delta_phi_fits
import mpetPorts

eps = -1e-12

# Define some variable types
mole_frac_t = daeVariableType(name="mole_frac_t", units=unit(),
        lowerBound=0, upperBound=1, initialGuess=0.25,
        absTolerance=1e-6)
elec_pot_t = daeVariableType(name="elec_pot_t", units=unit(),
        lowerBound=-1e20, upperBound=1e20, initialGuess=0,
        absTolerance=1e-5)

class mod0D1var(daeModel):
    def __init__(self, Name, Parent=None, Description="", ndD=None):
        daeModel.__init__(self, Name, Parent, Description)

class mod1D1var(daeModel):
    def __init__(self, Name, Parent=None, Description="", ndD=None):
        daeModel.__init__(self, Name, Parent, Description)

        if (ndD is None):
            raise Exception("Need input parameter dictionary")
        self.ndD = ndD

        # Domain
        self.Dmn = daeDomain("discretizationDomain", self, unit(),
                "discretization domain")

        # Variables
        self.c =  daeVariable("c", mole_frac_t, self,
                "Concentration in active particle",
                [self.Dmn])
        self.cbar = daeVariable("cbar", mole_frac_t, self,
                "Average concentration in active particle")
        self.dcbardt = daeVariable("dcbardt", no_t, self,
                "Rate of particle filling")
        if self.ndD["simSurfCond"]:
            self.phi = daeVariable("phi", elec_pot_t, self,
                    "Electric potential within the particle",
                    [self.Dmn])
        else:
            self.phi = False

        # Ports
#        self.portOut = mpetPorts.portFromParticle()
        self.portInLyte = mpetPorts.portFromElyte()
        self.portInBulk = mpetPorts.portFromBulk()
        self.phi_lyte = self.portInLyte.phi_lyte()
        self.c_lyte = self.portInLyte.c_lyte()
        self.phi_m = self.portInBulk.phi_m()

    def DeclareEquations(self):
        ndD = self.ndD
        N = ndD["N"]
        r_vec, volfrac_vec = get_unit_solid_discr(
                ndD['particleShape'],
                ndD['particleType'], N)

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

        # Equations for concentration evolution
        (Mmat, RHS_c) = self.calc_dcdt()
        dcdt_vec = np
        dcdt_vec = np.empty(Nij, dtype=object)
        dcdt_vec[0:N] = [self.c.dt(k) for k in range(Nij)]
        LHS_vec = MX(Mmat, dcdt_vec)
        for k in range(N):
            eq = self.CreateEquation("dcsdt_discr{k}".format(k=k))
            eq.Residual = LHS_vec[k] - RHS_c[k]
#            eq.Residual = (LHS_vec[k] - RHS_c[k] + noisevec[k]()) # NOISE

        # Equations for potential drop along particle, if desired
        if ndD['simSurfCond']:
            # Conservation of charge in the solid particles with
            # Ohm's Law
#            LHS = self.calc_part_surf_LHS()
            phi_tmp = np.empty(N + 2, dtype=object)
            phi_tmp[1:-1] = [self.phi(k) for k in range(N)]
            # BC's -- touching e- supply on both sides
            phi_tmp[0] = self.phi_m
            phi_tmp[-1] = self.phi_m
            dx = 1./N
            phi_edges (phi_tmp[0:-1] + phi_tmp[1:])/2.
            scond_vec = ndD["scond"] * np.exp(-1*(phi_edges - self.phi_m))
            curr_dens = -scond_vec * (np.diff(phi_tmp, 1) / dx)
            LHS = np.diff(curr_dens, 1)/dx
            k0_part = ndD["k0"]
            for k in range(N):
                eq = self.CreateEquation(
                        "charge_cons_discr{k}".format(
                            i=i,j=j,k=k,l=l))
                RHS = self.c.dt(k) / k0_part
                eq.Residual = LHS[k] - RHS

        for eq in self.Equations:
            eq.CheckUnitsConsistency = False

    def calc_sld_dcs_dt(self):
        # Get some useful information
        ndD = self.ndD
#        simSurfCond = ndD['simSurfCond']
        solidType = ndD['particleType']
        solidShape = ndD['particleShape']
#        rxnType = ndD['rxnType']
        delPhiEqFit = ndD['delPhiEqFit']
        # Get variables for this particle/electrode volume
        phi_lyte = self.phi_lyte
        phi_m = self.phi_m
        c_lyte = self.c_lyte
        cbar = self.cbar() # only used for ACR/CHR
        # Get the relevant parameters for this particle
        k0 = ndD["k0"]
        kappa = ndD["kappa"] # only used for ACR/CHR
        lmbda = ndD["lambda"] # Only used for Marcus/MHC
        alpha = ndD["alpha"] # Only used for BV
        Omga = ndD["Omga"]
        Ds = ndD["Dsld"] # Only used for "diffn" or "CHR"
        # We need the (non-dimensional) temperature to get the
        # reaction rate dependence correct
        T = ndD["T"]
        # Number of volumes in current particle
        N = ndD["N"]
        # Concentration (profile?) in the solid
        c = np.empty(N, dtype=object)
        c[:] = [self.c(k) for k in range(N)]
        # Assume dilute electrolyte
        act_O = c_lyte
        # If we're also simulating potential drop along the solid,
        # use that instead of self.phi_c(i)
        if ndD["simSurfCond"]:
            phi = np.empty(N, dtype=object)
            phi[:] = [self.phi(k) for k in range(N)]
        else:
            phi = phi_m
        mu_O = T*np.log(Max(eps, act_O)) + phi_lyte - phi

        # Calculate chemical potential of reduced state
#        if solidType in ["ACR", "homog", "homog_sdn"]:
        if solidType in ["ACR"]:
            # Make a blank array to allow for boundary conditions
            ctmp = np.empty(Nij+2, dtype=object)
            ctmp[1:-1] = c
            ctmp[0] = ndD["cwet"]
            ctmp[-1] = ndD["cwet"]
            dxs = 1./N
            curv = np.diff(ctmp, 2)/(dxs**2)
            mu_R = ( mu_reg_sln(c, Omga, T) - kappa*curv
                    + ndD["B"][l]*(c - cbar) )
            # XXX -- Temp dependence!
            act_R = np.exp(mu_R/T)
#            # eta = electrochem pot_R - electrochem pot_O
#            # eta = mu_R - mu_O
#            if delPhiEqFit:
#                material = ndD['material'][l]
#                fits = delta_phi_fits.DPhiFits(ndD["T"])
#                phifunc = fits.materialData[material]
#                delta_phi_eq = (phifunc(c[-1], ndD["dphi_eq_ref"])
#                        + np.log(act_O))
#                eta = (phi - phi_lyte) - delta_phi_eq
#            else:
#                eta = mu_R - mu_O
            eta = self.get_eta(c, act_O, mu_R, mu_O, T,
                    ndD["delPhiEqFit"], ndD["dphi_eq_ref"],
                    ndD["material"])
            Rxn = self.get_rxn_rate(eta, c, act_R, c_lyte, act_O, k0,
                    T, rxnType, lmbda, alpha)
            M = sprs.eye(N, N, format="csr")
            return (M, Rxn)

        elif solidType in ["diffn", "CHR"]:
            # For discretization background, see Zeng & Bazant 2013
            # Mass matrix is common for spherical shape, diffn or CHR
            Rs = 1. # (non-dimensionalized by itself)
            r_vec, volfrac_vec = get_unit_solid_discr(
                    solidShape, solidType, N)
            edges = np.hstack((0, (r_vec[0:-1] + r_vec[1:])/2, Rs))
            if solidShape == "sphere":
                Vp = 4./3. * np.pi * Rs**3
            elif solidShape == "cylinder":
                Vp = np.pi * Rs**2  # per unit height
            vol_vec = Vp * volfrac_vec
            dr = r_vec[1] - r_vec[0]
            M1 = sprs.diags([1./8, 3./4, 1./8], [-1, 0, 1],
                    shape=(N, N), format="csr")
            M1[1, 0] = M1[-2, -1] = 1./4
            M2 = sprs.diags(vol_vec, 0, format="csr")
            if solidShape == "sphere":
                M = M1*M2
            elif solidShape == "cylinder":
                M = M2

            # Chemical potentials, reactions first
            # Diffn
            if solidType in ["diffn"]:
#                Flux_vec = calc_Flux_diffn(c, Ds, Flux_bc, dr)
#                # Only sphere is implemented for diffusion right now
#                RHS = np.empty(Nij, dtype=object)
#                c_diffs = np.diff(c)
#                RHS[1:Nij - 1] = 4*np.pi*(
#                        Ds*(r_vec[1:Nij - 1] + dr/2.)**2*c_diffs[1:]/dr -
#                        Ds*(r_vec[1:Nij - 1] - dr/2.)**2*c_diffs[:-1]/dr )
#                RHS[0] = 4*np.pi*Ds*(dr/2.)**2*c_diffs[0]/dr
                # Take the surface concentration
                c_surf = c[-1]
                # Assuming dilute solid
                act_R_surf = c_surf
                mu_R_surf = T*np.log(Max(eps, act_R_surf))

            # CHR
            elif solidType in ["CHR"]:
                # mu_R is like for ACR, except kappa*curv term
                mu_R = ( mu_reg_sln(c, Omga, T) +
                        ndD["B"]*(c - cbar) )
                # Surface conc gradient given by natural BC
                beta_s = ndD["beta_s"]
                curv_c = calc_curv_c(c, dr, r_vec, Rs, beta_s,
                        particleShape)
                mu_R -= kappa*curv_c
#                if solidShape == "sphere":
#                    mu_R[0] -= 3 * kappa * (2*c[1] - 2*c[0])/dr**2
#                    mu_R[1:Nij - 1] -= kappa * (np.diff(c, 2)/dr**2 +
#                            (c[2:] - c[0:-2])/(dr*r_vec[1:-1])
#                            )
#                    mu_R[Nij - 1] -= kappa * ((2./Rs)*beta_s +
#                            (2*c[-2] - 2*c[-1] + 2*dr*beta_s)/dr**2
#                            )
#                elif solidShape == "cylinder":
#                    mu_R[0] -= 2 * kappa * (2*c[1] - 2*c[0])/dr**2
#                    mu_R[1:Nij - 1] -= kappa * (np.diff(c, 2)/dr**2 +
#                            (c[2:] - c[0:-2])/(2 * dr*r_vec[1:-1])
#                            )
#                    mu_R[Nij - 1] -= kappa * ((1./Rs)*beta_s +
#                            (2*c[-2] - 2*c[-1] + 2*dr*beta_s)/dr**2
#                            )
                mu_R_surf = mu_R[-1]
                # With chem potentials, can now calculate fluxes
#                Flux_vec = np.empty(Nij+1, dtype=object)
#                Flux_vec[0] = 0  # Symmetry at r=0
#                c_edges = (c_sld[0:-1]  + c_sld[1:])/2.
#                Flux_vec[1:Nij] = (Ds/T * (1 - c_edges) * c_edges *
#                        np.diff(mu_R)/dr)
                # Take the surface concentration
                c_surf = c_sld[-1]
                act_R_surf = np.exp(mu_R_surf/T)
            # Figure out overpotential
            eta = self.get_eta(c_surf, act_O, mu_R_surf, mu_O, T,
                    ndD["delPhiEqFit"], ndD["dphi_eq_ref"],
                    ndD["material"])
#            Rxn = self.get_rxn_rate(eta, c, act_R_surf, c_lyte, act_O)
            Rxn = self.get_rxn_rate(eta, c_surf, act_R_surf, c_lyte,
                    act_O, k0, T, rxnType, lmbda, alpha)
#            if delPhiEqFit:
#                material = ndD['material'][l]
#                fits = delta_phi_fits.DPhiFits(ndD["T"])
#                phifunc = fits.materialData[material]
#                delta_phi_eq = phifunc(c_surf, ndD["dphi_eq_ref"][l])
#                eta = (phi_m - phi_lyte) - delta_phi_eq
#            else:
#                eta = mu_R_surf - mu_O
#            # Calculate reaction rate
#            if rxnType == "Marcus":
#                Rxn = self.R_Marcus(k0, lmbda, c_lyte, c_surf, eta, T)
#            elif rxnType == "BV":
#                Rxn = self.R_BV(k0, alpha, c_surf, act_O, act_R_surf, eta, T)
#            elif rxnType == "MHC":
#                k0_MHC = k0/self.MHC_kfunc(0., lmbda)
#                Rxn = self.R_MHC(k0_MHC, lmbda, eta, T, c_surf, c_lyte)
            # Finish up RHS discretization at particle surface
            Flux_bc = ndD["delta_L"] * Rxn
            if solidType in ["diffn"]:
                Flux_vec = calc_Flux_diffn(c, Ds, Flux_bc, dr, T)
#                RHS[-1] = 4*np.pi*(Rs**2 * ndD["delta_L"][l][i, j] * Rxn -
#                        Ds*(Rs - dr/2)**2*c_diffs[-1]/dr )
            elif solidType in ["CHR"]:
#                Flux_vec[Nij] = ndD["delta_L"][l][i, j] * Rxn
                Flux_vec = calc_Flux_CHR(c, mu_R, Ds, Flux_bc, dr, T)
            if solidShape == "sphere":
                area_vec = 4*np.pi*edges**2
            elif solidShape == "cylinder":
                area_vec = 2*np.pi*edges  # per unit height
            RHS = np.diff(Flux_vec*area_vec)
            return (M, RHS)

def get_rxn_rate(eta, c_sld, act_R, c_lyte, act_O, k0, T, rxnType,
        lmbda=None, alpha=None):
#    ndD = self.ndD
#    rxnType = ndD["rxnType"]
#    k0 = ndD["k0"]
#    T = ndD["T"]
    if rxnType == "Marcus":
        Rate = R_Marcus(k0, lmbda, c_lyte, c_sld, eta, T)
    elif rxnType == "BV":
        Rate = R_BV(k0, alpha, c_sld, act_O, act_R, eta, T)
    elif rxnType == "MHC":
        k0_MHC = k0/MHC_kfunc(0., lmbda)
        Rate = R_MHC(k0_MHC, lmbda, eta, T, c_sld, c_lyte)
    return Rate

def get_eta(c, act_O, mu_R, mu_O, T, delPhiEqFit, dphi_eq_ref=None,
        material=None):
    if delPhiEqFit:
#        material = ndD['material'][l]
#        fits = delta_phi_fits.DPhiFits(ndD["T"])
#        delta_phi_eq = (phifunc(c, ndD["dphi_eq_ref"])
#                + np.log(act_O))
        fits = delta_phi_fits.DPhiFits(T)
        phifunc = fits.materialData[material]
        delta_phi_eq = (phifunc(c, dphi_eq_ref)
                + np.log(act_O))
        eta = (phi - phi_lyte) - delta_phi_eq
    else:
        eta = mu_R - mu_O
    return eta

def get_unit_solid_discr(solidShape, solidType, N):
    if solidShape == "C3" and solidType in ["ACR"]:
        r_vec = None
        # For 1D particle, the vol fracs are simply related to the
        # length discretization
        volfrac_vec = (1./N) * np.ones(N)  # scaled to 1D particle volume
        return r_vec, volfrac_vec
#    if solidType in ["homog", "homog_sdn"]:
#        r_vec = None
#        volfrac_vec = np.ones(1)
#        return r_vec, volfrac_vec
    if solidShape == "sphere":
        Rs = 1.
        dr = Rs/(N - 1)
        r_vec = np.linspace(0, Rs, N)
        vol_vec = 4*np.pi*(r_vec**2 * dr + (1./12)*dr**3)
        vol_vec[0] = 4*np.pi*(1./24)*dr**3
        vol_vec[-1] = (4./3)*np.pi*(Rs**3 - (Rs - dr/2.)**3)
        Vp = 4./3.*np.pi*Rs**3
        volfrac_vec = vol_vec/Vp
        return r_vec, volfrac_vec
    if solidShape == "cylinder":
        Rs = 1.
        h = 1.
        dr = Rs / (N - 1)
        r_vec = np.linspace(0, Rs, N)
        vol_vec = np.pi * h * 2 * r_vec * dr
        vol_vec[0] = np.pi * h * dr**2 / 4.
        vol_vec[-1] = np.pi * h * (Rs * dr - dr**2 / 4.)
        Vp = np.pi * Rs**2 * h
        volfrac_vec = vol_vec / Vp
        return r_vec, volfrac_vec
    else:
        raise NotImplementedError("Fix shape volumes!")

#def calc_Flux_diffn(c, Ds, Flux_bc, r_vec, dr):
def calc_Flux_diffn(c, Ds, Flux_bc, dr, T):
    N = len(c)
    Flux_vec = np.empty(N+1, dtype=object)
    Flux_vec[0] = 0 # Symmetry at r=0
    Flux_vec[-1] = Flux_bc
    Flux_vec[1:N] = Ds/T * np.diff(c)/dr
    return Flux_vec

def calc_Flux_CHR(c, mu, Ds, Flux_bc, dr, T):
    N = len(c)
    Flux_vec = np.empty(N+1, dtype=ojbect)
    Flux_vec[0] = 0 # Symmetry at r=0
    Flux_vec[-1] = Flux_bc
    c_edges = (c[0:-1] + c[1:])/2.
    # Keep the concentration between 0 and 1
    c_edges = np.array([Max(1e-6, c_edges[i]) for i in range(N)])
    c_edges = np.array([Min(1-1e-6, c_edges[i]) for i in range(N)])
    Flux_vec[1:N] = (Ds/T * (1-c_edges) * c_edges *
            np.diff(mu)/dr)
    return Flux_vec

def calc_curv_c(c, dr, r_vec, Rs, beta_s, particleShape):
    N = len(c)
    curv = np.empty(N, dtype=ojbect)
    if particleShape == "sphere":
        curv[0] = 3 * (2*c[1] - 2*c[0]) / dr**2
        curv[1:N-1] = (np.diff(c, 2)/dr**2 +
                (c[2:] - c[0:-2])/(dr*r_vec[1:-1]))
        curv[N-1] = ((2./Rs)*beta_s +
                (2*c[-2] - 2*c[-1] + 2*dr*beta_s)/dr**2)
    elif particleShape == "cylinder":
        curv[0] = 2 * (2*c[1] - 2*c[0]) / dr**2
        curv[1:N-1] = (np.diff(c, 2)/dr**2 +
                (c[2:] - c[0:-2])/(2 * dr*r_vec[1:-1]))
        curv[N-1] = ((1./Rs)*beta_s +
                (2*c[-2] - 2*c[-1] + 2*dr*beta_s)/dr**2)
    else:
        raise NotImplementedError("calc_curv_c only for sphere and cylinder")
    return curv

def mu_reg_sln(c, Omga, T):
    return np.array([ Omga*(1-2*c[i])
            + T*Log(Max(eps, c[i])/Max(eps, 1-c[i]))
            for i in range(len(c)) ])

def R_BV(k0, alpha, c_sld, act_O, act_R, eta, T):
    gamma_ts = (1./(1-c_sld))
    ecd = ( k0 * act_O**(1-alpha)
            * act_R**(alpha) / gamma_ts )
    Rate = ( ecd *
        (np.exp(-alpha*eta/T) - np.exp((1-alpha)*eta/T)) )
    return Rate

def R_Marcus(k0, lmbda, c_lyte, c_sld, eta, T):
    if type(c_sld) == np.ndarray:
        c_sld = np.array([Max(eps, c_sld[i]) for i in
            range(len(c_sld))])
    else:
        c_sld = Max(eps, c_sld)
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

def MHC_kfunc(eta, lmbda):
    a = 1. + np.sqrt(lmbda)
    if type(eta) == pyCore.adouble:
        ERF = Erf
    else:
        ERF = spcl.erf
    # evaluate with eta for oxidation, -eta for reduction
    return (np.sqrt(np.pi*lmbda) / (1 + np.exp(-eta))
            * (1. - ERF((lmbda - np.sqrt(a + eta**2))
                / (2*np.sqrt(lmbda)))))

def R_MHC(k0, lmbda, eta, T, c_sld, c_lyte):
    # See Zeng, Smith, Bai, Bazant 2014
    # Convert to "MHC overpotential"
    eta_f = eta + T*np.log(c_lyte/c_sld)
    gamma_ts = 1./(1. - c_sld)
    if type(eta) == np.ndarray:
        Rate = np.empty(len(eta), dtype=object)
        for i, etaval in enumerate(eta):
            krd = k0*MHC_kfunc(-eta_f[i], lmbda)
            kox = k0*MHC_kfunc(eta_f[i], lmbda)
            Rate[i] = (1./gamma_ts[i])*(krd*c_lyte - kox*c_sld[i])
    else:
        krd = k0*MHC_kfunc(-eta_f, lmbda)
        kox = k0*MHC_kfunc(eta_f, lmbda)
        Rate = (1./gamma_ts)*(krd*c_lyte - kox*c_sld)
    return Rate

def MX(mat, objvec):
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
