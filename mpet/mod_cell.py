"""The model defining the macroscopic cell.

This includes the equations defining
 - transport in the electrolyte
 - overall quantities such as current and voltage and their specification
 - overall filling fraction of each electrode
 - potential drop along the electrodes
 - potential drop between simulated particles
"""
import daetools.pyDAE as dae
from pyUnits import s

import numpy as np

import mpet.extern_funcs as extern_funcs
import mpet.geometry as geom
import mpet.mod_electrodes as mod_electrodes
import mpet.ports as ports
import mpet.props_elyte as props_elyte
import mpet.utils as utils
from mpet.config import constants
from mpet.daeVariableTypes import mole_frac_t, elec_pot_t, conc_t

# Dictionary of end conditions
endConditions = {
    1:"Vmax reached",
    2:"Vmin reached"}


class ModCell(dae.daeModel):
    def __init__(self, config, Name, Parent=None, Description=""):
        dae.daeModel.__init__(self, Name, Parent, Description)

        self.config = config
        self.profileType = config['profileType']
        Nvol = config["Nvol"]
        Npart = config["Npart"]
        self.trodes = trodes = config["trodes"]

        # Domains where variables are distributed
        self.DmnCell = {}  # domains over full cell dimensions
        self.DmnPart = {}  # domains over particles in each cell volume
        if config['have_separator']:  # If we have a separator
            self.DmnCell["s"] = dae.daeDomain(
                "DmnCell_s", self, dae.unit(),
                "Simulated volumes in the separator")
        for trode in trodes:
            self.DmnCell[trode] = dae.daeDomain(
                "DmnCell_{trode}".format(trode=trode), self, dae.unit(),
                "Simulated volumes in electrode {trode}".format(trode=trode))
            self.DmnPart[trode] = dae.daeDomain(
                "Npart_{trode}".format(trode=trode), self, dae.unit(),
                "Particles sampled in each control "
                + "volume in electrode {trode}".format(trode=trode))

        # Variables
        self.c_lyte = {}
        self.phi_lyte = {}
        self.phi_bulk = {}
        self.phi_part = {}
        self.R_Vp = {}
        self.ffrac = {}
        for trode in trodes:
            # Concentration/potential in electrode regions of elyte
            self.c_lyte[trode] = dae.daeVariable(
                "c_lyte_{trode}".format(trode=trode), conc_t, self,
                "Concentration in the elyte in electrode {trode}".format(trode=trode),
                [self.DmnCell[trode]])
            self.phi_lyte[trode] = dae.daeVariable(
                "phi_lyte_{trode}".format(trode=trode), elec_pot_t, self,
                "Electric potential in elyte in electrode {trode}".format(trode=trode),
                [self.DmnCell[trode]])
            self.phi_bulk[trode] = dae.daeVariable(
                "phi_bulk_{trode}".format(trode=trode), elec_pot_t, self,
                "Electrostatic potential in the bulk solid",
                [self.DmnCell[trode]])
            self.phi_part[trode] = dae.daeVariable(
                "phi_part_{trode}".format(trode=trode), elec_pot_t, self,
                "Electrostatic potential at each particle",
                [self.DmnCell[trode], self.DmnPart[trode]])
            self.R_Vp[trode] = dae.daeVariable(
                "R_Vp_{trode}".format(trode=trode), dae.no_t, self,
                "Rate of reaction of positives per electrode volume",
                [self.DmnCell[trode]])
            self.ffrac[trode] = dae.daeVariable(
                "ffrac_{trode}".format(trode=trode), mole_frac_t, self,
                "Overall filling fraction of solids in electrodes")
        if config['have_separator']:  # If we have a separator
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
        if ('a' not in config['trodes']) and (not config['have_separator']) and Nvol["c"] == 1:
            self.SVsim = True
        else:
            self.SVsim = False
        if not self.SVsim:
            # Ghost points (GP) to aid in boundary condition (BC) implemenation
            self.c_lyteGP_L = dae.daeVariable("c_lyteGP_L", conc_t, self, "c_lyte left BC GP")
            self.phi_lyteGP_L = dae.daeVariable(
                "phi_lyteGP_L", elec_pot_t, self, "phi_lyte left BC GP")
        self.phi_applied = dae.daeVariable(
            "phi_applied", elec_pot_t, self,
            "Overall battery voltage (at anode current collector)")
        self.phi_cell = dae.daeVariable(
            "phi_cell", elec_pot_t, self,
            "Voltage between electrodes (phi_applied less series resistance)")
        self.current = dae.daeVariable(
            "current", dae.no_t, self, "Total current of the cell")
        self.endCondition = dae.daeVariable(
            "endCondition", dae.no_t, self, "A nonzero value halts the simulation")

        # Create models for representative particles within electrode
        # volumes and ports with which to talk to them.
        self.portsOutLyte = {}
        self.portsOutBulk = {}
        self.particles = {}
        for trode in trodes:
            Nv = Nvol[trode]
            Np = Npart[trode]
            self.portsOutLyte[trode] = np.empty(Nv, dtype=object)
            self.portsOutBulk[trode] = np.empty((Nv, Np), dtype=object)
            self.particles[trode] = np.empty((Nv, Np), dtype=object)
            for vInd in range(Nv):
                self.portsOutLyte[trode][vInd] = ports.portFromElyte(
                    "portTrode{trode}vol{vInd}".format(trode=trode, vInd=vInd), dae.eOutletPort,
                    self, "Electrolyte port to particles")
                for pInd in range(Np):
                    self.portsOutBulk[trode][vInd,pInd] = ports.portFromBulk(
                        "portTrode{trode}vol{vInd}part{pInd}".format(
                            trode=trode, vInd=vInd, pInd=pInd),
                        dae.eOutletPort, self,
                        "Bulk electrode port to particles")
                    solidType = config[trode, "type"]
                    if solidType in constants.two_var_types:
                        pMod = mod_electrodes.Mod2var
                    elif solidType in constants.one_var_types:
                        pMod = mod_electrodes.Mod1var
                    else:
                        raise NotImplementedError("unknown solid type")
                    self.particles[trode][vInd,pInd] = pMod(
                        config, trode, vInd, pInd,
                        Name="partTrode{trode}vol{vInd}part{pInd}".format(
                            trode=trode, vInd=vInd, pInd=pInd),
                        Parent=self)
                    self.ConnectPorts(self.portsOutLyte[trode][vInd],
                                      self.particles[trode][vInd,pInd].portInLyte)
                    self.ConnectPorts(self.portsOutBulk[trode][vInd,pInd],
                                      self.particles[trode][vInd,pInd].portInBulk)

    def DeclareEquations(self):
        dae.daeModel.DeclareEquations(self)

        # Some values of domain lengths
        trodes = self.trodes
        config = self.config
        Nvol = config["Nvol"]
        Npart = config["Npart"]
        Nlyte = np.sum(list(Nvol.values()))

        # Define the overall filling fraction in the electrodes
        for trode in trodes:
            eq = self.CreateEquation("ffrac_{trode}".format(trode=trode))
            eq.Residual = self.ffrac[trode]()
            dx = 1./Nvol[trode]
            # Make a float of Vtot, total particle volume in electrode
            # Note: for some reason, even when "factored out", it's a bit
            # slower to use Sum(self.psd_vol_ac[l].array([], [])
            tmp = 0
            for vInd in range(Nvol[trode]):
                for pInd in range(Npart[trode]):
                    Vj = config["psd_vol_FracVol"][trode][vInd,pInd]
                    tmp += self.particles[trode][vInd,pInd].cbar() * Vj * dx
            eq.Residual -= tmp

        # Define dimensionless R_Vp for each electrode volume
        for trode in trodes:
            for vInd in range(Nvol[trode]):
                eq = self.CreateEquation(
                    "R_Vp_trode{trode}vol{vInd}".format(vInd=vInd, trode=trode))
                # Start with no reaction, then add reactions for each
                # particle in the volume.
                RHS = 0
                # sum over particle volumes in given electrode volume
                for pInd in range(Npart[trode]):
                    # The volume of this particular particle
                    Vj = config["psd_vol_FracVol"][trode][vInd,pInd]
                    RHS += -(config["beta"][trode] * (1-config["poros"][trode])
                             * config["P_L"][trode] * Vj
                             * self.particles[trode][vInd,pInd].dcbardt())
                eq.Residual = self.R_Vp[trode](vInd) - RHS

        # Define output port variables
        for trode in trodes:
            for vInd in range(Nvol[trode]):
                eq = self.CreateEquation(
                    "portout_c_trode{trode}vol{vInd}".format(vInd=vInd, trode=trode))
                eq.Residual = (self.c_lyte[trode](vInd)
                               - self.portsOutLyte[trode][vInd].c_lyte())
                eq = self.CreateEquation(
                    "portout_p_trode{trode}vol{vInd}".format(vInd=vInd, trode=trode))
                phi_lyte = self.phi_lyte[trode](vInd)
                eq.Residual = (phi_lyte - self.portsOutLyte[trode][vInd].phi_lyte())
                for pInd in range(Npart[trode]):
                    eq = self.CreateEquation(
                        "portout_pm_trode{trode}v{vInd}p{pInd}".format(
                            vInd=vInd, pInd=pInd, trode=trode))
                    eq.Residual = (self.phi_part[trode](vInd, pInd)
                                   - self.portsOutBulk[trode][vInd,pInd].phi_m())

            # Simulate the potential drop along the bulk electrode
            # solid phase
            simBulkCond = config['simBulkCond'][trode]
            if simBulkCond:
                # Calculate the RHS for electrode conductivity
                phi_tmp = utils.add_gp_to_vec(utils.get_var_vec(self.phi_bulk[trode], Nvol[trode]))
                porosvec = utils.pad_vec(utils.get_const_vec(
                    (1-self.config["poros"][trode])**(1-config["BruggExp"][trode]), Nvol[trode]))
                poros_walls = utils.mean_harmonic(porosvec)
                if trode == "a":  # anode
                    # Potential at the current collector is from
                    # simulation
                    phi_tmp[0] = self.phi_cell()
                    # No current passes into the electrolyte
                    phi_tmp[-1] = phi_tmp[-2]
                else:  # cathode
                    phi_tmp[0] = phi_tmp[1]
                    # Potential at current at current collector is
                    # reference (set)
                    phi_tmp[-1] = config["phi_cathode"]
                dx = config["L"][trode]/Nvol[trode]
                dvg_curr_dens = np.diff(-poros_walls*config["sigma_s"][trode]
                                        * np.diff(phi_tmp)/dx)/dx
            # Actually set up the equations for bulk solid phi
            for vInd in range(Nvol[trode]):
                eq = self.CreateEquation(
                    "phi_ac_trode{trode}vol{vInd}".format(vInd=vInd, trode=trode))
                if simBulkCond:
                    eq.Residual = -dvg_curr_dens[vInd] - self.R_Vp[trode](vInd)
                else:
                    if trode == "a":  # anode
                        eq.Residual = self.phi_bulk[trode](vInd) - self.phi_cell()
                    else:  # cathode
                        eq.Residual = self.phi_bulk[trode](vInd) - config["phi_cathode"]

            # Simulate the potential drop along the connected
            # particles
            simPartCond = config['simPartCond'][trode]
            for vInd in range(Nvol[trode]):
                phi_bulk = self.phi_bulk[trode](vInd)
                for pInd in range(Npart[trode]):
                    G_l = config["G"][trode][vInd,pInd]
                    phi_n = self.phi_part[trode](vInd, pInd)
                    if pInd == 0:  # reference bulk phi
                        phi_l = phi_bulk
                    else:
                        phi_l = self.phi_part[trode](vInd, pInd-1)
                    if pInd == (Npart[trode] - 1):  # No particle at end of "chain"
                        G_r = 0
                        phi_r = phi_n
                    else:
                        G_r = config["G"][trode][vInd,pInd+1]
                        phi_r = self.phi_part[trode](vInd, pInd+1)
                    # charge conservation equation around this particle
                    eq = self.CreateEquation(
                        "phi_ac_trode{trode}vol{vInd}part{pInd}".format(
                            vInd=vInd, trode=trode, pInd=pInd))
                    if simPartCond:
                        # -dcsbar/dt = I_l - I_r
                        eq.Residual = (
                            self.particles[trode][vInd,pInd].dcbardt()
                            + ((-G_l * (phi_n - phi_l))
                               - (-G_r * (phi_r - phi_n))))
                    else:
                        eq.Residual = self.phi_part[trode](vInd, pInd) - phi_bulk

        # If we have a single electrode volume (in a perfect bath),
        # electrolyte equations are simple
        if self.SVsim:
            eq = self.CreateEquation("c_lyte")
            eq.Residual = self.c_lyte["c"].dt(0) - 0
            eq = self.CreateEquation("phi_lyte")
            eq.Residual = self.phi_lyte["c"](0) - self.phi_cell()
        else:
            disc = geom.get_elyte_disc(Nvol, config["L"], config["poros"], config["BruggExp"])
            cvec = utils.get_asc_vec(self.c_lyte, Nvol)
            dcdtvec = utils.get_asc_vec(self.c_lyte, Nvol, dt=True)
            phivec = utils.get_asc_vec(self.phi_lyte, Nvol)
            Rvvec = utils.get_asc_vec(self.R_Vp, Nvol)
            # Apply concentration and potential boundary conditions
            # Ghost points on the left and no-gradients on the right
            ctmp = np.hstack((self.c_lyteGP_L(), cvec, cvec[-1]))
            phitmp = np.hstack((self.phi_lyteGP_L(), phivec, phivec[-1]))

            Nm_edges, i_edges = get_lyte_internal_fluxes(ctmp, phitmp, disc, config)

            # If we don't have a porous anode:
            # 1) the total current flowing into the electrolyte is set
            # 2) assume we have a Li foil with BV kinetics and the specified rate constant
            eqC = self.CreateEquation("GhostPointC_L")
            eqP = self.CreateEquation("GhostPointP_L")
            if 'a' not in config["trodes"]:
                # Concentration BC from mass flux
                eqC.Residual = Nm_edges[0]
                # Phi BC from BV at the foil
                # We assume BV kinetics with alpha = 0.5,
                # exchange current density, ecd = k0_foil * c_lyte**(0.5)
                cWall = .5*(ctmp[0] + ctmp[1])
                ecd = config["k0_foil"]*cWall**0.5
                # -current = ecd*(exp(-eta/2) - exp(eta/2))
                # note negative current because positive current is
                # oxidation here
                # -current = ecd*(-2*sinh(eta/2))
                # eta = 2*arcsinh(-current/(-2*ecd))
                BVfunc = -self.current() / ecd
                eta_eff = 2*np.arcsinh(-BVfunc/2.)
                eta = eta_eff + self.current()*config["Rfilm_foil"]
#                # Infinitely fast anode kinetics
#                eta = 0.
                # eta = mu_R - mu_O = -mu_O (evaluated at interface)
                # mu_O = [T*ln(c) +] phiWall - phi_cell = -eta
                # phiWall = -eta + phi_cell [- T*ln(c)]
                phiWall = -eta + self.phi_cell()
                if config["elyteModelType"] == "dilute":
                    phiWall -= config["T"]*np.log(cWall)
                # phiWall = 0.5 * (phitmp[0] + phitmp[1])
                eqP.Residual = phiWall - .5*(phitmp[0] + phitmp[1])
            # We have a porous anode -- no flux of charge or anions through current collector
            else:
                eqC.Residual = ctmp[0] - ctmp[1]
                eqP.Residual = phitmp[0] - phitmp[1]

            dvgNm = np.diff(Nm_edges)/disc["dxvec"]
            dvgi = np.diff(i_edges)/disc["dxvec"]
            for vInd in range(Nlyte):
                # Mass Conservation (done with the anion, although "c" is neutral salt conc)
                eq = self.CreateEquation("lyte_mass_cons_vol{vInd}".format(vInd=vInd))
                eq.Residual = disc["porosvec"][vInd]*dcdtvec[vInd] + (1./config["num"])*dvgNm[vInd]
                # Charge Conservation
                eq = self.CreateEquation("lyte_charge_cons_vol{vInd}".format(vInd=vInd))
                eq.Residual = -dvgi[vInd] + config["zp"]*Rvvec[vInd]

        # Define the total current. This must be done at the capacity
        # limiting electrode because currents are specified in
        # C-rates.
        eq = self.CreateEquation("Total_Current")
        eq.Residual = self.current()
        limtrode = config["limtrode"]
        dx = 1./Nvol[limtrode]
        rxn_scl = config["beta"][limtrode] * (1-config["poros"][limtrode]) \
            * config["P_L"][limtrode]
        for vInd in range(Nvol[limtrode]):
            if limtrode == "a":
                eq.Residual -= dx * self.R_Vp[limtrode](vInd)/rxn_scl
            else:
                eq.Residual += dx * self.R_Vp[limtrode](vInd)/rxn_scl
        # Define the measured voltage, offset by the "applied" voltage
        # by any series resistance.
        # phi_cell = phi_applied - I*R
        eq = self.CreateEquation("Measured_Voltage")
        eq.Residual = self.phi_cell() - (
            self.phi_applied() - config["Rser"]*self.current())

        if self.profileType == "CC":
            # Total Current Constraint Equation
            eq = self.CreateEquation("Total_Current_Constraint")
            if config["tramp"] > 0:
                eq.Residual = self.current() - (
                    config["currPrev"] + (config["currset"] - config["currPrev"])
                    * (1 - np.exp(-dae.Time()/(config["tend"]*config["tramp"]))))
            else:
                eq.Residual = self.current() - config["currset"]
        elif self.profileType == "CV":
            # Keep applied potential constant
            eq = self.CreateEquation("applied_potential")
            if config["tramp"] > 0:
                eq.Residual = self.phi_applied() - (
                    config["phiPrev"] + (config["Vset"] - config["phiPrev"])
                    * (1 - np.exp(-dae.Time()/(config["tend"]*config["tramp"])))
                    )
            else:
                eq.Residual = self.phi_applied() - config["Vset"]
        elif self.profileType == "CP":
            # constant power constraint
            ndDVref = config["c", "phiRef"]
            if 'a' in config["trodes"]:
                ndDVref = config["c", "phiRef"] - config["a", "phiRef"]
            eq = self.CreateEquation("Total_Power_Constraint")
            # adding Vref since P = V*I
            if config["tramp"] > 0:
                eq.Residual = self.current()*(self.phi_applied() + ndDVref) - (
                    config["currPrev"]*(config["phiPrev"] + ndDVref)
                    + (config["power"] - (config["currPrev"]*(config["phiPrev"]
                                                              + ndDVref)))
                    * (1 - np.exp(-dae.Time()/(config["tend"]*config["tramp"])))
                    )
            else:
                eq.Residual = self.current()*(self.phi_applied() + ndDVref) - config["power"]
        elif self.profileType == "CCsegments":
            if config["tramp"] > 0:
                config["segments_setvec"][0] = config["currPrev"]
                self.segSet = extern_funcs.InterpTimeScalar(
                    "segSet", self, dae.unit(), dae.Time(),
                    config["segments_tvec"], config["segments_setvec"])
                eq = self.CreateEquation("Total_Current_Constraint")
                eq.Residual = self.current() - self.segSet()

            # CCsegments implemented as discontinuous equations
            else:
                # First segment
                time = config["segments"][0][1]
                self.IF(dae.Time() < dae.Constant(time*s), 1.e-3)
                eq = self.CreateEquation("Total_Current_Constraint")
                eq.Residual = self.current() - config["segments"][0][0]

                # Middle segments
                for i in range(1,len(config["segments"])-1):
                    time = time+config["segments"][i][1]
                    self.ELSE_IF(dae.Time() < dae.Constant(time*s), 1.e-3)
                    eq = self.CreateEquation("Total_Current_Constraint")
                    eq.Residual = self.current() - config["segments"][i][0]

                # Last segment
                self.ELSE()
                eq = self.CreateEquation("Total_Current_Constraint")
                eq.Residual = self.current() - config["segments"][-1][0]
                self.END_IF()

        elif self.profileType == "CVsegments":
            if config["tramp"] > 0:
                config["segments_setvec"][0] = config["phiPrev"]
                self.segSet = extern_funcs.InterpTimeScalar(
                    "segSet", self, dae.unit(), dae.Time(),
                    config["segments_tvec"], config["segments_setvec"])
                eq = self.CreateEquation("applied_potential")
                eq.Residual = self.phi_applied() - self.segSet()

            # CVsegments implemented as discontinuous equations
            else:
                # First segment
                time = config["segments"][0][1]
                self.IF(dae.Time() < dae.Constant(time*s), 1.e-3)
                eq = self.CreateEquation("applied_potential")
                eq.Residual = self.phi_applied() - config["segments"][0][0]

                # Middle segments
                for i in range(1,len(config["segments"])-1):
                    time = time+config["segments"][i][1]
                    self.ELSE_IF(dae.Time() < dae.Constant(time*s), 1.e-3)
                    eq = self.CreateEquation("applied_potential")
                    eq.Residual = self.phi_applied() - config["segments"][i][0]

                # Last segment
                self.ELSE()
                eq = self.CreateEquation("applied_potential")
                eq.Residual = self.phi_applied() - config["segments"][-1][0]
                self.END_IF()

        for eq in self.Equations:
            eq.CheckUnitsConsistency = False

        # Ending conditions for the simulation
        if self.profileType in ["CC", "CCsegments"]:
            # Vmax reached
            self.ON_CONDITION((self.phi_applied() <= config["phimin"])
                              & (self.endCondition() < 1),
                              setVariableValues=[(self.endCondition, 1)])

            # Vmin reached
            self.ON_CONDITION((self.phi_applied() >= config["phimax"])
                              & (self.endCondition() < 1),
                              setVariableValues=[(self.endCondition, 2)])


def get_lyte_internal_fluxes(c_lyte, phi_lyte, disc, config):
    zp, zm, nup, num = config["zp"], config["zm"], config["nup"], config["num"]
    nu = nup + num
    T = config["T"]
    dxd1 = disc["dxd1"]
    eps_o_tau = disc["eps_o_tau"]

    # Get concentration at cell edges using weighted mean
    wt = utils.pad_vec(disc["dxvec"])
    c_edges_int = utils.weighted_linear_mean(c_lyte, wt)

    if config["elyteModelType"] == "dilute":
        # Get porosity at cell edges using weighted harmonic mean
        eps_o_tau_edges = utils.weighted_linear_mean(eps_o_tau, wt)
        Dp = eps_o_tau_edges * config["Dp"]
        Dm = eps_o_tau_edges * config["Dm"]
#        Np_edges_int = nup*(-Dp*np.diff(c_lyte)/dxd1
#                            - Dp*zp*c_edges_int*np.diff(phi_lyte)/dxd1)
        Nm_edges_int = num*(-Dm*np.diff(c_lyte)/dxd1
                            - Dm/T*zm*c_edges_int*np.diff(phi_lyte)/dxd1)
        i_edges_int = (-((nup*zp*Dp + num*zm*Dm)*np.diff(c_lyte)/dxd1)
                       - (nup*zp**2*Dp + num*zm**2*Dm)/T*c_edges_int*np.diff(phi_lyte)/dxd1)
#        i_edges_int = zp*Np_edges_int + zm*Nm_edges_int
    elif config["elyteModelType"] == "SM":
        D_fs, sigma_fs, thermFac, tp0 = getattr(props_elyte,config["SMset"])()[:-1]

        # Get diffusivity and conductivity at cell edges using weighted harmonic mean
        D_edges = utils.weighted_harmonic_mean(eps_o_tau*D_fs(c_lyte, T), wt)
        sigma_edges = utils.weighted_harmonic_mean(eps_o_tau*sigma_fs(c_lyte, T), wt)

        sp, n = config["sp"], config["n"]
        # there is an error in the MPET paper, temperature dependence should be
        # in sigma and not outside of sigma
        i_edges_int = -sigma_edges * (
            np.diff(phi_lyte)/dxd1
            + nu*T*(sp/(n*nup)+tp0(c_edges_int, T)/(zp*nup))
            * thermFac(c_edges_int, T)
            * np.diff(np.log(c_lyte))/dxd1
            )
        Nm_edges_int = num*(-D_edges*np.diff(c_lyte)/dxd1
                            + (1./(num*zm)*(1-tp0(c_edges_int, T))*i_edges_int))
    return Nm_edges_int, i_edges_int
