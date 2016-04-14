import numpy as np
import scipy.sparse as sprs
import scipy.special as spcl

import daetools.pyDAE as dae

import muRfuncs
import mpetPorts
import externFuncs

eps = -1e-12

class mod2var(dae.daeModel):
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

        # Define some variable types
        atol = ndD_s["absTol"]
        mole_frac_t = dae.daeVariableType(
                name="mole_frac_t", units=dae.unit(), lowerBound=0,
                upperBound=1, initialGuess=0.25, absTolerance=atol)
        elec_pot_t = dae.daeVariableType(
                name="elec_pot_t", units=dae.unit(), lowerBound=-1e20,
                upperBound=1e20, initialGuess=0, absTolerance=atol)
        # Variables
        self.c1 =  dae.daeVariable("c1", mole_frac_t, self,
                "Concentration in 'layer' 1 of active particle",
                [self.Dmn])
        self.c2 =  dae.daeVariable("c2", mole_frac_t, self,
                "Concentration in 'layer' 2 of active particle",
                [self.Dmn])
        self.cbar = dae.daeVariable("cbar", mole_frac_t, self,
                "Average concentration in active particle")
        self.c1bar = dae.daeVariable("c1bar", mole_frac_t, self,
                "Average concentration in 'layer' 1 of active particle")
        self.c2bar = dae.daeVariable("c2bar", mole_frac_t, self,
                "Average concentration in 'layer' 2 of active particle")
        self.dcbardt = dae.daeVariable("dcbardt", dae.no_t, self,
                "Rate of particle filling")

        # Ports
        self.portInLyte = mpetPorts.portFromElyte("portInLyte",
                dae.eInletPort, self, "Inlet port from electrolyte")
        self.portInBulk = mpetPorts.portFromBulk("portInBulk",
                dae.eInletPort, self, "Inlet port from e- conducting phase")
        self.phi_lyte = self.portInLyte.phi_lyte()
        self.c_lyte = self.portInLyte.c_lyte()
        self.phi_m = self.portInBulk.phi_m()

    def DeclareEquations(self):
        dae.daeModel.DeclareEquations(self)
        ndD = self.ndD
        N = ndD["N"] # number of grid points in particle
        T = self.ndD_s["T"] # nondimensional temperature
        r_vec, volfrac_vec = get_unit_solid_discr(ndD['shape'], N)

        # Prepare the Ideal Solution log ratio terms
        self.ISfuncs1 = self.ISfuncs2 = None
        if ndD["logPad"]:
            self.ISfuncs1 = np.array(
                    [externFuncs.LogRatio("LR1", self, dae.unit(),
                        self.c1(k)) for k in range(N)])
            self.ISfuncs2 = np.array(
                    [externFuncs.LogRatio("LR2", self, dae.unit(),
                        self.c2(k)) for k in range(N)])
        ISfuncs = (self.ISfuncs1, self.ISfuncs2)

        # Prepare noise
        self.noise1 = self.noise2 = None
        if ndD["noise"]:
            numnoise = 200
            noise_prefac = 1e-7
            tvec = np.linspace(0., 1.05*self.ndD_s["tend"], numnoise)
            noise_data1 = noise_prefac*np.random.randn(numnoise, N)
            noise_data2 = noise_prefac*np.random.randn(numnoise, N)
            # Previous_output is common for all external functions
            previous_output1 = []
            previous_output2 = []
            self.noise1 = [externFuncs.InterpTimeVector("noise1",
                self, dae.unit(), dae.Time(), tvec, noise_data1,
                previous_output1, _position_) for _position_ in
                range(N)]
            self.noise2 = [externFuncs.InterpTimeVector("noise2",
                self, dae.unit(), dae.Time(), tvec, noise_data2,
                previous_output2, _position_) for _position_ in
                range(N)]
        noises = (self.noise1, self.noise2)

        # Figure out mu_O, mu of the oxidized state
        mu_O, act_lyte = calc_mu_O(self.c_lyte, self.phi_lyte, self.phi_m, T,
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
        eq.Residual = self.cbar() - 0.5*(self.c1bar() + self.c2bar())

        # Define average rate of filling of particle
        eq = self.CreateEquation("dcbardt")
        eq.Residual = self.dcbardt()
        for k in range(N):
            eq.Residual -= 0.5*(self.c1.dt(k) + self.c2.dt(k)) * volfrac_vec[k]

        c1 = np.empty(N, dtype=object)
        c2 = np.empty(N, dtype=object)
        c1[:] = [self.c1(k) for k in range(N)]
        c2[:] = [self.c2(k) for k in range(N)]
        if ndD["type"] in ["diffn2", "CHR2"]:
            # Equations for 1D particles of 1 field varible
            self.sldDynamics1D2var(c1, c2, mu_O, act_lyte, ISfuncs, noises)
        elif ndD["type"] in ["homog2", "homog2_sdn"]:
            # Equations for 0D particles of 1 field variables
            self.sldDynamics0D2var(c1, c2, mu_O, act_lyte, ISfuncs, noises)

        for eq in self.Equations:
            eq.CheckUnitsConsistency = False

    def sldDynamics0D2var(self, c1, c2, muO, act_lyte, ISfuncs, noises):
        ndD = self.ndD
        N = ndD["N"]
        T = self.ndD_s["T"]
        c1_surf = c1
        c2_surf = c2
        (mu1R_surf, mu2R_surf), (act1R_surf, act2R_surf) = (
                calc_muR(
                    (c1_surf, c2_surf),
                    (self.c1bar(), self.c2bar()),
                    T, ndD, ISfuncs))
        eta1 = calc_eta(mu1R_surf, muO)
        eta2 = calc_eta(mu2R_surf, muO)
        Rxn1 = calc_rxn_rate(eta1, c1_surf, self.c_lyte, ndD["k0"],
                T, ndD["rxnType"], act1R_surf, act_lyte, ndD["lambda"],
                ndD["alpha"])
        Rxn2 = calc_rxn_rate(eta2, c2_surf, self.c_lyte, ndD["k0"],
                T, ndD["rxnType"], act2R_surf, act_lyte, ndD["lambda"],
                ndD["alpha"])

        dc1dt_vec = np.empty(N, dtype=object)
        dc2dt_vec = np.empty(N, dtype=object)
        dc1dt_vec[0:N] = [self.c1.dt(k) for k in range(N)]
        dc2dt_vec[0:N] = [self.c2.dt(k) for k in range(N)]
        LHS1_vec = dc1dt_vec
        LHS2_vec = dc2dt_vec
        noise1, noise2 = noises
        for k in range(N):
            eq1 = self.CreateEquation("dc1sdt")
            eq2 = self.CreateEquation("dc2sdt")
            eq1.Residual = LHS1_vec[k] - Rxn1[k]
            eq2.Residual = LHS2_vec[k] - Rxn2[k]
            if ndD["noise"]:
                eq1.Residual += noise1[k]()
                eq2.Residual += noise2[k]()
        return

    def sldDynamics1D2var(self, c1, c2, muO, act_lyte, ISfuncs, noises):
        ndD = self.ndD
        N = ndD["N"]
        T = self.ndD_s["T"]
        # Equations for concentration evolution
        # Mass matrix, M, where M*dcdt = RHS, where c and RHS are vectors
        Mmat = get_Mmat(ndD['shape'], N)
        dr, edges = get_dr_edges(ndD['shape'], N)

        # Get solid particle chemical potential, overpotential, reaction rate
        c1_surf = mu1_R_surf = act1_R_surf = None
        c2_surf = mu2_R_surf = act2_R_surf = None
        if ndD["type"] in ["diffn2", "CHR2"]:
            (mu1R, mu2R), (act1R, act2R) = calc_muR((c1, c2),
                    (self.c1bar(), self.c2bar()), T, ndD, ISfuncs)
            c1_surf = c1[-1]
            c2_surf = c2[-1]
            mu1R_surf, act1R_surf = mu1R[-1], act1R[-1]
            mu2R_surf, act2R_surf = mu2R[-1], act2R[-1]
        eta1 = calc_eta(mu1R_surf, muO)
        eta2 = calc_eta(mu2R_surf, muO)
        Rxn1 = calc_rxn_rate(eta1, c1_surf, self.c_lyte, ndD["k0"],
                T, ndD["rxnType"], act1R_surf, act_lyte, ndD["lambda"],
                ndD["alpha"])
        Rxn2 = calc_rxn_rate(eta2, c2_surf, self.c_lyte, ndD["k0"],
                T, ndD["rxnType"], act2R_surf, act_lyte, ndD["lambda"],
                ndD["alpha"])

        # Get solid particle fluxes (if any) and RHS
        if ndD["type"] in ["diffn2", "CHR2"]:
            # Positive reaction (reduction, intercalation) is negative
            # flux of Li at the surface.
            Flux1_bc = -0.5 * ndD["delta_L"] * Rxn1
            Flux2_bc = -0.5 * ndD["delta_L"] * Rxn2
            if ndD["type"] == "diffn2":
                Flux1_vec, Flux2_vec = calc_Flux_diffn2(c1, c2,
                        ndD["D"], Flux1_bc, Flux2_bc, dr, T)
            elif ndD["type"] == "CHR2":
                Flux1_vec, Flux2_vec = calc_Flux_CHR2(c1, c2, mu1R, mu2R,
                        ndD["D"], Flux1_bc, Flux2_bc, dr, T)
            if ndD["shape"] == "sphere":
                area_vec = 4*np.pi*edges**2
            elif ndD["shape"] == "cylinder":
                area_vec = 2*np.pi*edges  # per unit height
            RHS1 = -np.diff(Flux1_vec * area_vec)
            RHS2 = -np.diff(Flux2_vec * area_vec)
#            kinterlayer = 1e2
#            interLayerRxn = (kinterlayer * (1 - c1_sld) *
#                    (1 - c2_sld) * (act1_R - act2_R))
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
        return

class mod1var(dae.daeModel):
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

        # Define some variable types
        atol = ndD_s["absTol"]
        mole_frac_t = dae.daeVariableType(
                name="mole_frac_t", units=dae.unit(), lowerBound=0,
                upperBound=1, initialGuess=0.25, absTolerance=atol)
        elec_pot_t = dae.daeVariableType(
                name="elec_pot_t", units=dae.unit(), lowerBound=-1e20,
                upperBound=1e20, initialGuess=0, absTolerance=atol)
        # Variables
        self.c =  dae.daeVariable("c", mole_frac_t, self,
                "Concentration in active particle",
                [self.Dmn])
        self.cbar = dae.daeVariable("cbar", mole_frac_t, self,
                "Average concentration in active particle")
        self.dcbardt = dae.daeVariable("dcbardt", dae.no_t, self,
                "Rate of particle filling")

        # Ports
        self.portInLyte = mpetPorts.portFromElyte("portInLyte",
                dae.eInletPort, self, "Inlet port from electrolyte")
        self.portInBulk = mpetPorts.portFromBulk("portInBulk",
                dae.eInletPort, self, "Inlet port from e- conducting phase")
        self.phi_lyte = self.portInLyte.phi_lyte()
        self.c_lyte = self.portInLyte.c_lyte()
        self.phi_m = self.portInBulk.phi_m()

    def DeclareEquations(self):
        dae.daeModel.DeclareEquations(self)
        ndD = self.ndD
        N = ndD["N"] # number of grid points in particle
        T = self.ndD_s["T"] # nondimensional temperature
        r_vec, volfrac_vec = get_unit_solid_discr(ndD['shape'], N)

        # Prepare the Ideal Solution log ratio terms
        self.ISfuncs = None
        if ndD["logPad"]:
            self.ISfuncs = np.array(
                    [externFuncs.LogRatio("LR", self, dae.unit(),
                        self.c(k)) for k in range(N)])

        # Prepare noise
        self.noise = None
        if ndD["noise"]:
            numnoise = 200
            noise_prefac = 1e-5
            tvec = np.linspace(0., 1.05*self.ndD_s["tend"], numnoise)
            noise_data = noise_prefac*np.random.randn(numnoise, N)
            # Previous_output is common for all external functions
            previous_output = []
            self.noise = [externFuncs.InterpTimeVector("noise", self,
                dae.unit(), dae.Time(), tvec, noise_data, previous_output,
                _position_) for _position_ in range(N)]

        # Figure out mu_O, mu of the oxidized state
        mu_O, act_lyte = calc_mu_O(self.c_lyte, self.phi_lyte, self.phi_m, T,
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
            self.sldDynamics1D1var(c, mu_O, act_lyte, self.ISfuncs,
                    self.noise)
        elif ndD["type"] in ["homog", "homog_sdn"]:
            # Equations for 0D particles of 1 field variables
            self.sldDynamics0D1var(c, mu_O, act_lyte, self.ISfuncs,
                    self.noise)

        for eq in self.Equations:
            eq.CheckUnitsConsistency = False

    def sldDynamics0D1var(self, c, muO, act_lyte, ISfuncs, noise):
        ndD = self.ndD
        N = ndD["N"]
        T = self.ndD_s["T"]
        c_surf = c
        muR_surf, actR_surf = calc_muR(c_surf, self.cbar(), T, ndD, ISfuncs)
        eta = calc_eta(muR_surf, muO)
        Rxn = calc_rxn_rate(eta, c_surf, self.c_lyte, ndD["k0"],
                T, ndD["rxnType"], actR_surf, act_lyte, ndD["lambda"],
                ndD["alpha"])

        dcdt_vec = np.empty(N, dtype=object)
        dcdt_vec[0:N] = [self.c.dt(k) for k in range(N)]
        LHS_vec = dcdt_vec
        for k in range(N):
            eq = self.CreateEquation("dcsdt")
            eq.Residual = LHS_vec[k] - Rxn[k]
            if ndD["noise"]:
                eq.Residual += noise[k]()
        return

    def sldDynamics1D1var(self, c, muO, act_lyte, ISfuncs, noise):
        ndD = self.ndD
        N = ndD["N"]
        T = self.ndD_s["T"]
        # Equations for concentration evolution
        # Mass matrix, M, where M*dcdt = RHS, where c and RHS are vectors
        Mmat = get_Mmat(ndD['shape'], N)
        dr, edges = get_dr_edges(ndD['shape'], N)

        # Get solid particle chemical potential, overpotential, reaction rate
        c_surf = mu_R_surf = act_R_surf = None
        if ndD["type"] in ["ACR"]:
            c_surf = c
            muR_surf, actR_surf = calc_muR(c_surf, self.cbar(), T, ndD, ISfuncs)
        elif ndD["type"] in ["diffn", "CHR"]:
            muR, actR = calc_muR(c, self.cbar(), T, ndD, ISfuncs)
            c_surf = c[-1]
            muR_surf = muR[-1]
            actR_surf = actR[-1]
        eta = calc_eta(muR_surf, muO)
        Rxn = calc_rxn_rate(eta, c_surf, self.c_lyte, ndD["k0"],
                T, ndD["rxnType"], actR_surf, act_lyte, ndD["lambda"],
                ndD["alpha"])

        # Get solid particle fluxes (if any) and RHS
        if ndD["type"] in ["ACR"]:
            RHS = Rxn
        elif ndD["type"] in ["diffn", "CHR"]:
            # Positive reaction (reduction, intercalation) is negative
            # flux of Li at the surface.
            Flux_bc = -ndD["delta_L"] * Rxn
            if ndD["type"] == "diffn":
                Flux_vec = calc_Flux_diffn(c, ndD["D"], Flux_bc, dr, T)
            elif ndD["type"] == "CHR":
                Flux_vec = calc_Flux_CHR(c, muR, ndD["D"], Flux_bc, dr, T)
            if ndD["shape"] == "sphere":
                area_vec = 4*np.pi*edges**2
            elif ndD["shape"] == "cylinder":
                area_vec = 2*np.pi*edges  # per unit height
            RHS = -np.diff(Flux_vec * area_vec)

        dcdt_vec = np.empty(N, dtype=object)
        dcdt_vec[0:N] = [self.c.dt(k) for k in range(N)]
        LHS_vec = MX(Mmat, dcdt_vec)
        for k in range(N):
            eq = self.CreateEquation("dcsdt_discr{k}".format(k=k))
            eq.Residual = LHS_vec[k] - RHS[k]
            if ndD["noise"]:
                eq.Residual += noise[k]()

        return

def calc_rxn_rate(eta, c_sld, c_lyte, k0, T, rxnType,
        act_R=None, act_lyte=None, lmbda=None, alpha=None):
    if rxnType == "Marcus":
        Rate = R_Marcus(k0, lmbda, c_lyte, c_sld, eta, T)
    elif rxnType[0:2] == "BV":
        Rate = R_BV(k0, alpha, c_lyte, c_sld, act_lyte, act_R,
                eta, T, rxnType)
    elif rxnType == "MHC":
        k0_MHC = k0/MHC_kfunc(0., lmbda)
        Rate = R_MHC(k0_MHC, lmbda, eta, T, c_sld, c_lyte)
    return Rate

def calc_eta(muR, muO):
    return muR - muO

def get_unit_solid_discr(Shape, N):
    if N == 1: # homog particle, hopefully
        r_vec = None
        volfrac_vec = np.ones(1)
    elif Shape == "C3":
        r_vec = None
        # For 1D particle, the vol fracs are simply related to the
        # length discretization
        volfrac_vec = (1./N) * np.ones(N)  # scaled to 1D particle volume
    elif Shape == "sphere":
        Rs = 1. # (non-dimensionalized by itself)
        dr = Rs/(N - 1)
        r_vec = np.linspace(0, Rs, N)
        vol_vec = 4*np.pi*(r_vec**2 * dr + (1./12)*dr**3)
        vol_vec[0] = 4*np.pi*(1./24)*dr**3
        vol_vec[-1] = (4./3)*np.pi*(Rs**3 - (Rs - dr/2.)**3)
        Vp = 4./3.*np.pi*Rs**3
        volfrac_vec = vol_vec/Vp
    elif Shape == "cylinder":
        Rs = 1. # (non-dimensionalized by itself)
        h = 1.
        dr = Rs / (N - 1)
        r_vec = np.linspace(0, Rs, N)
        vol_vec = np.pi * h * 2 * r_vec * dr
        vol_vec[0] = np.pi * h * dr**2 / 4.
        vol_vec[-1] = np.pi * h * (Rs * dr - dr**2 / 4.)
        Vp = np.pi * Rs**2 * h
        volfrac_vec = vol_vec / Vp
    else:
        raise NotImplementedError("Fix shape volumes!")
    return r_vec, volfrac_vec

def get_dr_edges(shape, N):
    r_vec = get_unit_solid_discr(shape, N)[0]
    dr = edges = None
    if r_vec is not None:
        Rs = 1.
        dr = r_vec[1] - r_vec[0]
        edges = np.hstack((0, (r_vec[0:-1] + r_vec[1:])/2, Rs))
    return dr, edges

def get_Mmat(shape, N):
    r_vec, volfrac_vec = get_unit_solid_discr(shape, N)
    if shape == "C3":
        Mmat = sprs.eye(N, N, format="csr")
    elif shape in ["sphere", "cylinder"]:
        dr = r_vec[1] - r_vec[0]
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
        M1[1, 0] = M1[-2, -1] = 1./4
        M2 = sprs.diags(vol_vec, 0, format="csr")
        Mmat = M1*M2
    return Mmat

def calc_Flux_diffn(c, Ds, Flux_bc, dr, T):
    N = len(c)
    Flux_vec = np.empty(N+1, dtype=object)
    Flux_vec[0] = 0 # Symmetry at r=0
    Flux_vec[-1] = Flux_bc
    Flux_vec[1:N] = -Ds/T * np.diff(c)/dr
    return Flux_vec

def calc_Flux_CHR(c, mu, Ds, Flux_bc, dr, T):
    N = len(c)
    Flux_vec = np.empty(N+1, dtype=object)
    Flux_vec[0] = 0 # Symmetry at r=0
    Flux_vec[-1] = Flux_bc
    c_edges = 2*(c[0:-1] * c[1:])/(c[0:-1] + c[1:])
    # Keep the concentration between 0 and 1
    c_edges = np.array([Max(1e-6, c_edges[i]) for i in range(N-1)])
    c_edges = np.array([Min(1-1e-6, c_edges[i]) for i in range(N-1)])
    Flux_vec[1:N] = -(Ds/T * (1-c_edges) * c_edges *
            np.diff(mu)/dr)
    return Flux_vec

def calc_Flux_CHR2(c1, c2, mu1_R, mu2_R, Ds, Flux1_bc, Flux2_bc, dr, T):
    if isinstance(c1[0], dae.pyCore.adouble):
        MIN, MAX = Min, Max
    else:
        MIN, MAX = min, max
    N = len(c1)
    Flux1_vec = np.empty(N+1, dtype=object)
    Flux2_vec = np.empty(N+1, dtype=object)
    Flux1_vec[0] = 0. # symmetry at r=0
    Flux2_vec[0] = 0. # symmetry at r=0
    Flux1_vec[-1] = Flux1_bc
    Flux2_vec[-1] = Flux2_bc
    c1_edges = 2*(c1[0:-1] * c1[1:])/(c1[0:-1] + c1[1:])
    c2_edges = 2*(c2[0:-1] * c2[1:])/(c2[0:-1] + c2[1:])
    # keep the concentrations between 0 and 1
    c1_edges = np.array([MAX(1e-6, c1_edges[i]) for i in
            range(len(c1_edges))])
    c1_edges = np.array([MIN((1-1e-6), c1_edges[i]) for i in
            range(len(c1_edges))])
    c2_edges = np.array([MAX(1e-6, c2_edges[i]) for i in
            range(len(c1_edges))])
    c2_edges = np.array([MIN((1-1e-6), c2_edges[i]) for i in
            range(len(c1_edges))])
    cbar_edges = 0.5*(c1_edges + c2_edges)
    Flux1_vec[1:N] = -(Ds/T * (1 - c1_edges)**(1.0) * c1_edges *
            np.diff(mu1_R)/dr)
    Flux2_vec[1:N] = -(Ds/T * (1 - c2_edges)**(1.0) * c2_edges *
            np.diff(mu2_R)/dr)
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
    muRfunc = muRfuncs.muRfuncs(T, ndD).muRfunc
    muR_ref = ndD["muR_ref"]
    muR, actR = muRfunc(c, cbar, muR_ref, ISfuncs)
    return muR, actR

def R_BV(k0, alpha, c_lyte, c_sld, act_lyte, act_R, eta, T, rxnType):
    if rxnType in ["BV", "BV_gMod01"]:
        if act_R is None:
            act_R = c_sld/(1-c_sld)
        if rxnType == "BV":
            gamma_ts = (1./(1-c_sld))
        elif rxnType == "BV_gMod01":
            gamma_ts = (1./(c_sld*(1-c_sld)))
        ecd = ( k0 * act_lyte**(1-alpha)
                * act_R**(alpha) / gamma_ts )
    elif rxnType in ["BV_raw", "BV_mod01", "BV_mod02"]:
        if rxnType == "BV_raw":
            ecd = k0
        elif rxnType == "BV_mod01":
            ecd = ( k0 * c_lyte**(1-alpha)
                    * (1.0 - c_sld)**(1 - alpha) * c_sld**alpha )
        elif rxnType == "BV_mod02":
            ecd = ( k0 * c_lyte**(1-alpha)
                    * (0.5 - c_sld)**(1 - alpha) * c_sld**alpha )
    Rate = ecd * (np.exp(-alpha*eta/T) - np.exp((1-alpha)*eta/T))
    return Rate

def R_Marcus(k0, lmbda, c_lyte, c_sld, eta, T):
    if isinstance(c_sld, np.ndarray):
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
    if isinstance(eta, dae.pyCore.adouble):
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
    if isinstance(eta, np.ndarray):
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
                    mat.data[low:up] * objvec[mat.indices[low:up]] )
        else:
            out[i] = 0.0
    return out
