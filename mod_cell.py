import daetools.pyDAE as dae
import numpy as np

import extern_funcs
import mod_electrodes
import ports
import props_elyte


class ModCell(dae.daeModel):
    def __init__(self, Name, Parent=None, Description="", ndD_s=None,
                 ndD_e=None):
        dae.daeModel.__init__(self, Name, Parent, Description)

        if (ndD_s is None) or (ndD_e is None):
            raise Exception("Need input parameter dictionaries")
        self.ndD = ndD_s
        self.profileType = ndD_s['profileType']
        Nvol = ndD_s["Nvol"]
        Npart = ndD_s["Npart"]
        self.trodes = trodes = ndD_s["trodes"]

        # Domains where variables are distributed
        self.DmnCell = {}  # domains over full cell dimensions
        self.DmnPart = {}  # domains over particles in each cell volume
        if Nvol["s"] >= 1:  # If we have a separator
            self.DmnCell["s"] = dae.daeDomain(
                "DmnCell_s", self, dae.unit(),
                "Simulated volumes in the separator")
        for l in trodes:
            self.DmnCell[l] = dae.daeDomain(
                "DmnCell_{l}".format(l=l), self, dae.unit(),
                "Simulated volumes in electrode {l}".format(l=l))
            self.DmnPart[l] = dae.daeDomain(
                "Npart_{l}".format(l=l), self, dae.unit(),
                "Particles sampled in each control " +
                "volume in electrode {l}".format(l=l))
            Nv = Nvol[l]
            Np = Npart[l]

        # Define some variable types
        atol = ndD_s["absTol"]
        mole_frac_t = dae.daeVariableType(
            name="mole_frac_t", units=dae.unit(), lowerBound=0,
            upperBound=1, initialGuess=0.25, absTolerance=atol)
        conc_t = dae.daeVariableType(
            name="conc_t", units=dae.unit(), lowerBound=0,
            upperBound=1e20, initialGuess=1.00, absTolerance=atol)
        elec_pot_t = dae.daeVariableType(
            name="elec_pot_t", units=dae.unit(), lowerBound=-1e20,
            upperBound=1e20, initialGuess=0, absTolerance=atol)
        # Variables
        self.c_lyte = {}
        self.phi_lyte = {}
        self.phi_bulk = {}
        self.phi_part = {}
        self.j_plus = {}
        self.ffrac = {}
        for l in trodes:
            # Concentration/potential in electrode regions of elyte
            self.c_lyte[l] = dae.daeVariable(
                "c_lyte_{l}".format(l=l), conc_t, self,
                "Concentration in the elyte in electrode {l}".format(l=l),
                [self.DmnCell[l]])
            self.phi_lyte[l] = dae.daeVariable(
                "phi_lyte_{l}".format(l=l), elec_pot_t, self,
                "Electric potential in elyte in electrode {l}".format(l=l),
                [self.DmnCell[l]])
            self.phi_bulk[l] = dae.daeVariable(
                "phi_bulk_{l}".format(l=l), elec_pot_t, self,
                "Electrostatic potential in the bulk solid",
                [self.DmnCell[l]])
            self.phi_part[l] = dae.daeVariable(
                "phi_part_{l}".format(l=l), elec_pot_t, self,
                "Electrostatic potential at each particle",
                [self.DmnCell[l], self.DmnPart[l]])
            self.j_plus[l] = dae.daeVariable(
                "j_plus_{l}".format(l=l), dae.no_t, self,
                "Rate of reaction of positives per solid volume",
                [self.DmnCell[l]])
            self.ffrac[l] = dae.daeVariable(
                "ffrac_{l}".format(l=l), mole_frac_t, self,
                "Overall filling fraction of solids in electrodes")
        if Nvol["s"] >= 1:  # If we have a separator
            self.c_lyte["s"] = dae.daeVariable(
                "c_lyte_s", conc_t, self,
                "Concentration in the electrolyte in the separator",
                [self.DmnCell["s"]])
            self.phi_lyte["s"] = dae.daeVariable(
                "phi_lyte_s", elec_pot_t, self,
                "Electrostatic potential in electrolyte in separator",
                [self.DmnCell["s"]])
        # Note if we're doing a single electrode volume simulation
        # It will be in a perfect bath of electrolyte at the applied
        # potential.
        if Nvol["a"] == 0 and Nvol["s"] == 0 and Nvol["c"] == 1:
            self.SVsim = True
        else:
            self.SVsim = False
        if not self.SVsim:
            self.c_lyteGP = dae.daeVariable(
                "c_lyteGP", conc_t, self,
                "Concentration in the electrolyte in " +
                "the boundary condition ghost point")
            self.phi_lyteGP = dae.daeVariable(
                "phi_lyteGP", elec_pot_t, self,
                "Electrostatic potential in electrolyte in " +
                "the boundary condition ghost point")
        self.phi_applied = dae.daeVariable(
            "phi_applied", elec_pot_t, self,
            "Overall battery voltage (at anode current collector)")
        self.phi_cell = dae.daeVariable(
            "phi_cell", elec_pot_t, self,
            "Voltage between electrodes (phi_applied less series resistance)")
        self.current = dae.daeVariable(
            "current", dae.no_t, self, "Total current of the cell")
        self.dummyVar = dae.daeVariable("dummyVar", dae.no_t, self, "dummyVar")

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
                self.portsOutLyte[l][i] = ports.portFromElyte(
                    "portTrode{l}vol{i}".format(l=l, i=i), dae.eOutletPort,
                    self, "Electrolyte port to particles")
                for j in range(Np):
                    self.portsOutBulk[l][i,j] = ports.portFromBulk(
                        "portTrode{l}vol{i}part{j}".format(l=l, i=i, j=j),
                        dae.eOutletPort, self,
                        "Bulk electrode port to particles")
                    solidType = ndD_e[l]["indvPart"][i,j]['type']
                    if solidType in ndD_s["2varTypes"]:
                        pMod = mod_electrodes.mod2var
                    elif solidType in ndD_s["1varTypes"]:
                        pMod = mod_electrodes.mod1var
                    else:
                        raise NotImplementedError("unknown solid type")
                    self.particles[l][i,j] = pMod(
                        "partTrode{l}vol{i}part{j}".format(l=l, i=i, j=j),
                        self, ndD=ndD_e[l]["indvPart"][i,j],
                        ndD_s=ndD_s)
                    self.ConnectPorts(self.portsOutLyte[l][i],
                                      self.particles[l][i,j].portInLyte)
                    self.ConnectPorts(self.portsOutBulk[l][i,j],
                                      self.particles[l][i,j].portInBulk)

    def DeclareEquations(self):
        dae.daeModel.DeclareEquations(self)

        # Some values of domain lengths
        trodes = self.trodes
        ndD = self.ndD
        Nvol = ndD["Nvol"]
        Npart = ndD["Npart"]
        Nlyte = np.sum(list(Nvol.values()))

        # Define the overall filling fraction in the electrodes
        for l in trodes:
            eq = self.CreateEquation("ffrac_{l}".format(l=l))
            eq.Residual = self.ffrac[l]()
            dx = 1./Nvol[l]
            # Make a float of Vtot, total particle volume in electrode
            # Note: for some reason, even when "factored out", it's a bit
            # slower to use Sum(self.psd_vol_ac[l].array([], [])
            tmp = 0
            for i in range(Nvol[l]):
                for j in range(Npart[l]):
                    Vj = ndD["psd_vol_FracVol"][l][i,j]
                    tmp += self.particles[l][i,j].cbar() * Vj * dx
            eq.Residual -= tmp

        # Define dimensionless j_plus for each electrode volume
        for l in trodes:
            for i in range(Nvol[l]):
                eq = self.CreateEquation(
                    "j_plus_trode{l}vol{i}".format(i=i, l=l))
                # Start with no reaction, then add reactions for each
                # particle in the volume.
                res = 0
                # sum over particle volumes in given electrode volume
                for j in range(Npart[l]):
                    # The volume of this particular particle
                    Vj = ndD["psd_vol_FracVol"][l][i,j]
                    res += self.particles[l][i,j].dcbardt() * Vj
                eq.Residual = self.j_plus[l](i) - res

        # Define output port variables
        for l in trodes:
            for i in range(Nvol[l]):
                eq = self.CreateEquation(
                    "portout_c_trode{l}vol{i}".format(i=i, l=l))
                eq.Residual = (self.c_lyte[l](i)
                               - self.portsOutLyte[l][i].c_lyte())
                eq = self.CreateEquation(
                    "portout_p_trode{l}vol{i}".format(i=i, l=l))
                phi_lyte = self.phi_lyte[l](i)
                eq.Residual = (phi_lyte - self.portsOutLyte[l][i].phi_lyte())
                for j in range(Npart[l]):
                    eq = self.CreateEquation(
                        "portout_pm_trode{l}v{i}p{j}".format(i=i, j=j, l=l))
                    eq.Residual = (self.phi_part[l](i, j) -
                                   self.portsOutBulk[l][i,j].phi_m())

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
                if l == "a":  # anode
                    # Potential at the current collector is from
                    # simulation
                    phi_tmp[0] = self.phi_cell()
                    # No current passes into the electrolyte
                    phi_tmp[-1] = phi_tmp[-2]
                else:  # cathode
                    phi_tmp[0] = phi_tmp[1]
                    # Potential at current at current collector is
                    # reference (set)
                    phi_tmp[-1] = ndD["phi_cathode"]
                dx = 1./Nvol[l]
                RHS_phi_tmp = -np.diff(
                    -porosvec*ndD["mcond"][l]*np.diff(phi_tmp)/dx)/dx
            # Actually set up the equations for bulk solid phi
            for i in range(Nvol[l]):
                eq = self.CreateEquation(
                    "phi_ac_trode{l}vol{i}".format(i=i, l=l))
                if simBulkCond:
                    eq.Residual = (-ndD["epsbeta"][l]*self.j_plus[l](i)
                                   - RHS_phi_tmp[i])
                else:
                    if l == "a":  # anode
                        eq.Residual = self.phi_bulk[l](i) - self.phi_cell()
                    else:  # cathode
                        eq.Residual = self.phi_bulk[l](i) - ndD["phi_cathode"]

            # Simulate the potential drop along the connected
            # particles
            simPartCond = ndD['simPartCond'][l]
            for i in range(Nvol[l]):
                phi_bulk = self.phi_bulk[l](i)
                for j in range(Npart[l]):
                    G_l = ndD["G"][l][i,j]
                    phi_n = self.phi_part[l](i, j)
                    if j == 0:  # reference bulk phi
                        phi_l = phi_bulk
                    else:
                        phi_l = self.phi_part[l](i, j-1)
                    if j == (Npart[l] - 1):  # No particle at end of "chain"
                        G_r = 0
                        phi_r = phi_n
                    else:
                        G_r = ndD["G"][l][i,j+1]
                        phi_r = self.phi_part[l](i, j+1)
                    # charge conservation equation around this particle
                    eq = self.CreateEquation(
                        "phi_ac_trode{l}vol{i}part{j}".format(i=i, l=l, j=j))
                    if simPartCond:
                        # -dcsbar/dt = I_l - I_r
                        eq.Residual = (
                            self.particles[l][i,j].dcbardt()
                            + ((-G_l * (phi_n - phi_l))
                               - (-G_r * (phi_r - phi_n))))
                    else:
                        eq.Residual = self.phi_part[l](i, j) - phi_bulk

        # If we have a single electrode volume (in a perfect bath),
        # electrolyte equations are simple
        if self.SVsim:
            eq = self.CreateEquation("c_lyte")
            eq.Residual = self.c_lyte["c"].dt(0) - 0
            eq = self.CreateEquation("phi_lyte")
            eq.Residual = self.phi_lyte["c"](0) - self.phi_cell()
        else:
            # Calculate RHS for electrolyte equations
            c_lyte = np.empty(Nlyte, dtype=object)
            c_lyte[0:Nvol["a"]] = [
                self.c_lyte["a"](i) for i in range(Nvol["a"])]  # anode
            c_lyte[Nvol["a"]:Nvol["a"] + Nvol["s"]] = [
                self.c_lyte["s"](i) for i in range(Nvol["s"])]  # separator
            c_lyte[Nvol["a"] + Nvol["s"]:Nlyte] = [
                self.c_lyte["c"](i) for i in range(Nvol["c"])]  # cathode
            phi_lyte = np.empty(Nlyte, dtype=object)
            phi_lyte[0:Nvol["a"]] = [
                self.phi_lyte["a"](i) for i in range(Nvol["a"])]  # anode
            phi_lyte[Nvol["a"]:Nvol["a"] + Nvol["s"]] = [
                self.phi_lyte["s"](i) for i in range(Nvol["s"])]  # separator
            phi_lyte[Nvol["a"] + Nvol["s"]:Nlyte] = [
                self.phi_lyte["c"](i) for i in range(Nvol["c"])]  # cathode
            (RHS_c, RHS_phi) = self.calc_lyte_RHS(
                c_lyte, phi_lyte, Nvol, Nlyte)
            # Equations governing the electrolyte in the separator
            offset = Nvol["a"]
            for i in range(Nvol["s"]):
                # Mass Conservation
                eq = self.CreateEquation(
                    "sep_lyte_mass_cons_vol{i}".format(i=i))
                eq.Residual = (ndD["poros"]["s"]*self.c_lyte["s"].dt(i)
                               - RHS_c[offset + i])
                # Charge Conservation
                eq = self.CreateEquation(
                    "sep_lyte_charge_cons_vol{i}".format(i=i))
                eq.Residual = (RHS_phi[offset + i])
            # Equations governing the electrolyte in the electrodes.
            # Here, we are coupled to the total reaction rates in the
            # solids.
            for l in trodes:
                if l == "a":  # anode
                    offset = 0
                else:  # cathode
                    offset = Nvol["a"] + Nvol["s"]
                for i in range(Nvol[l]):
                    # Mass Conservation
                    eq = self.CreateEquation(
                        "lyteMassCons_trode{l}vol{i}".format(i=i, l=l))
                    if ndD["elyteModelType"] == "dilute":
                        eq.Residual = (
                            ndD["poros"][l]*self.c_lyte[l].dt(i)
                            + ndD["epsbeta"][l]*(1-ndD["tp"])*self.j_plus[l](i)
                            - RHS_c[offset + i])
                    elif ndD["elyteModelType"] == "SM":
                        eq.Residual = (
                            ndD["poros"][l]*self.c_lyte[l].dt(i)
                            - RHS_c[offset + i])
                    # Charge Conservation
                    eq = self.CreateEquation(
                        "lyteChargeCons_trode{l}vol{i}".format(i=i, l=l))
                    eq.Residual = (ndD["epsbeta"][l]*self.j_plus[l](i)
                                   - RHS_phi[offset + i])

        # Define the total current. This must be done at the capacity
        # limiting electrode because currents are specified in
        # C-rates.
        eq = self.CreateEquation("Total_Current")
        eq.Residual = self.current()
        limtrode = ("c" if ndD["z"] < 1 else "a")
        dx = 1./Nvol[limtrode]
        for i in range(Nvol[limtrode]):
            if limtrode == "a":
                eq.Residual += dx * self.j_plus[limtrode](i)
            else:
                eq.Residual -= dx * self.j_plus[limtrode](i)
        # Define the measured voltage, offset by the "applied" voltage
        # by any series resistance.
        # phi_cell = phi_applied - I*R
        eq = self.CreateEquation("Measured_Voltage")
        eq.Residual = self.phi_cell() - (
            self.phi_applied() - ndD["Rser"]*self.current())

        if self.profileType == "CC":
            # Total Current Constraint Equation
            eq = self.CreateEquation("Total_Current_Constraint")
            eq.Residual = self.current() - (
                ndD["currPrev"] + (ndD["currset"] - ndD["currPrev"])
                * (1 - np.exp(-dae.Time()/(ndD["tend"]*ndD["tramp"])))
                )
        elif self.profileType == "CV":
            # Keep applied potential constant
            eq = self.CreateEquation("applied_potential")
            eq.Residual = self.phi_applied() - (
                ndD["phiPrev"] + (ndD["Vset"] - ndD["phiPrev"])
                * (1 - np.exp(-dae.Time()/(ndD["tend"]*ndD["tramp"])))
#                * 1
#                * np.tanh(dae.Time()/(45.0)))
                )
        elif "segments" in self.profileType:
            if self.profileType == "CCsegments":
                ndD["segments_setvec"][0] = ndD["currPrev"]
            elif self.profileType == "CVsegments":
                ndD["segments_setvec"][0] = ndD["phiPrev"]
            self.segSet = extern_funcs.InterpTimeScalar(
                "segSet", self, dae.unit(), dae.Time(),
                ndD["segments_tvec"], ndD["segments_setvec"])
            if self.profileType == "CCsegments":
                eq = self.CreateEquation("Total_Current_Constraint")
                eq.Residual = self.current() - self.segSet()
            elif self.profileType == "CVsegments":
                eq = self.CreateEquation("applied_potential")
                eq.Residual = self.phi_applied() - self.segSet()

        for eq in self.Equations:
            eq.CheckUnitsConsistency = False

        if self.profileType == "CC":
            # Set the condition to terminate the simulation upon reaching
            # a cutoff voltage.
            self.stopCondition = (
                ((self.phi_applied() <= ndD["phimin"])
                    | (self.phi_applied() >= ndD["phimax"]))
                & (self.dummyVar() < 1))
            self.ON_CONDITION(self.stopCondition,
                              setVariableValues=[(self.dummyVar, 2)])

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
        dxvec = np.empty(np.sum(list(Nvol.values())) + 2, dtype=object)
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
        porosvec[0:Nvol["a"]+1] = [
            ndD["poros"]["a"]**(3./2) for i in range(Nvol["a"]+1)]  # anode
        porosvec[Nvol["a"]+1:Nvol["a"]+1 + Nvol["s"]] = [
            ndD["poros"]["s"]**(3./2) for i in range(Nvol["s"])]  # separator
        porosvec[Nvol["a"]+1 + Nvol["s"]:] = [
            ndD["poros"]["c"]**(3./2) for i in range(Nvol["c"]+1)]  # cathode
        poros_edges = ((2*porosvec[1:]*porosvec[:-1])
                       / (porosvec[1:] + porosvec[:-1] + 1e-20))

        modType = ndD["elyteModelType"]
        if modType == "SM":
            D, kappa, thermFac, tp0 = props_elyte.getProps(ndD["SMset"])[:-1]
        # Apply concentration and potential boundary conditions
        ctmp = np.empty(Nlyte + 2, dtype=object)
        ctmp[1:-1] = cvec
        ctmp[0] = self.c_lyteGP()
        phitmp = np.empty(Nlyte + 2, dtype=object)
        phitmp[1:-1] = phivec
        phitmp[0] = self.phi_lyteGP()
        # If we don't have a real anode:
        # 1) the total current flowing into the electrolyte is set
        # 2) assume we have a Li foil with BV kinetics and the
        # specified rate constant
        eqC = self.CreateEquation("GhostPointC")
        eqP = self.CreateEquation("GhostPointP")

        if Nvol["a"] == 0:
            # Concentration BC from mass flux
            cWall = 2*ctmp[0]*ctmp[1]/(ctmp[0] + ctmp[1] + 1e-20)
            if modType == "dilute":
                Dwall, tp0wall = 1., ndD["tp"]
            elif modType == "SM":
                Dwall, tp0wall = D(cWall), tp0(cWall)
            currWall = self.current()*ndD["epsbeta"][limtrode]
            eqC.Residual = (poros_edges[0]*Dwall*(ctmp[1]-ctmp[0])/dxvec[0]
                            + currWall*(1-tp0wall))

            # Phi BC from BV at the foil
            # We assume BV kinetics with alpha = 0.5,
            # exchange current density, ecd = k0_foil * c_lyte**(0.5)
            ecd = ndD["k0_foil"]*cWall**0.5
            # -current = ecd*(exp(-eta/2) - exp(eta/2))
            # note negative current because positive current is
            # oxidation here
            # -current = ecd*(-2*sinh(eta/2))
            # eta = 2*arcsinh(-current/(-2*ecd))
            BVfunc = -self.current() / ecd
            eta_eff = 2*np.arcsinh(-BVfunc/2.)
            eta = eta_eff + self.current()*ndD["Rfilm_foil"]
#            # Infinitely fast anode kinetics
#            eta = 0.
            # eta = mu_R - mu_O = -mu_O (evaluated at interface)
            # mu_O = [T*ln(c) +] phiWall - phi_cell = -eta
            # phiWall = -eta + phi_cell [- T*ln(c)]
            phiWall = -eta + self.phi_cell()
            if modType == "dilute":
                phiWall -= ndD["T"]*np.log(cWall)
            # phiWall = 0.5 * (phitmp[0] + phitmp[1])
            eqP.Residual = phiWall - 0.5*(phitmp[0] + phitmp[1])
        else:  # porous anode -- no elyte flux at anode current collector
            eqC.Residual = ctmp[0] - ctmp[1]
            eqP.Residual = phitmp[0] - phitmp[1]
        # No electrolyte flux at the cathode current collector
        ctmp[-1] = ctmp[-2]
        phitmp[-1] = phitmp[-2]

        # We need average values of c_lyte for the current densities
        # at the finite volume boundaries
        c_edges = (2*ctmp[:-1]*ctmp[1:])/(ctmp[:-1] + ctmp[1:]+1e-20)
        # current density in the electrolyte
        if modType == "dilute":
            Dp = ndD["Dp"]
            Dm = ndD["Dm"]
            i_edges = poros_edges * (
                -((Dp - Dm)*np.diff(ctmp)/dxd1)
                - (zp*Dp - zm*Dm)*c_edges*np.diff(phitmp)/dxd1)
        elif modType == "SM":
            i_edges = -poros_edges * kappa(c_edges) * (
                np.diff(phitmp)/dxd1 -
                nu/nup*(1-tp0(c_edges)) *
                thermFac(c_edges) *
                np.diff(np.log(ctmp))/dxd1
                )
        # RHS for charge conservation equation
        RHS_phi = -np.diff(i_edges)/dxd2

        # RHS for species conservation equation
        if modType == "dilute":
            # Diffusive flux in the electrolyte
            cflux = -poros_edges*np.diff(ctmp)/dxd1
            # Divergence of the flux
            RHS_c = -np.diff(cflux)/dxd2
        elif modType == "SM":
            RHS_c = (
                np.diff(poros_edges * D(c_edges) *
                        np.diff(ctmp)/dxd1)/dxd2
                + 1./(zp*nup)*np.diff((1-tp0(c_edges))*i_edges)/dxd2
                )
        return (RHS_c, RHS_phi)
