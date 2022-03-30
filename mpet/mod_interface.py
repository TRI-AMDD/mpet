import numpy as np

import daetools.pyDAE as dae
from mpet import ports, utils, props_elyte
import mpet.geometry as geom
from mpet.daeVariableTypes import mole_frac_t, elec_pot_t


"""
Model for the interface between the electrolyte and active particles
"""


class InterfaceRegion(dae.daeModel):
    def __init__(self, config, Name, Parent=None, Description="",
                 cell=None, particle=None, vInd=None, pInd=None,
                 trode=None):
        super().__init__(Name, Parent, Description)
        self.config = config

        # Domain
        self.Dmn = dae.daeDomain("discretizationDomain", self, dae.unit(),
                                 "discretization domain")

        # Variables
        self.c = dae.daeVariable("c", mole_frac_t, self,
                                 "Concentration in interface",
                                 [self.Dmn])

        self.phi = dae.daeVariable("phi", elec_pot_t, self,
                                   "Electrical potential in interface",
                                   [self.Dmn])

        # Ports
        self.portInLyte = ports.portFromElyte(
            "portInLyte", dae.eInletPort, self,
            "Inlet port from electrolyte")

        self.portInParticle = ports.portFromParticle(
            "portInParticle", dae.eInletPort, self,
            "Inlet port from particle")

        # Note: using portFromElyte here because the port from the
        # interface to the particles transfers the same variables as
        # the port from the elyte to the interface, hence the same
        # class can be used
        self.portOutInterfaceParticle = ports.portFromElyte(
            "portOutInterfaceParticle", dae.eOutletPort, self,
            "Port from interface to particles")

        self.portOutInterfaceElyte = ports.portFromInterface(
            "portOutInterfaceElyte", dae.eOutletPort, self,
            "Port from interface to elyte")

        # Particle
        self.particle = particle

        # Cell
        self.cell = cell

        # Volume and particle indices
        self.vInd = vInd
        self.pInd = pInd
        self.trode = trode

    def DeclareEquations(self):
        super().DeclareEquations()
        config = self.config
        Nvol = config["Nvol_i"]

        disc = geom.get_interface_disc(Nvol, config["L_i"],
                                       config["poros_i"], config["BruggExp_i"])
        cvec = utils.get_var_vec(self.c, Nvol)
        dcdtvec = utils.get_var_vec(self.c, Nvol, dt=True)
        phivec = utils.get_var_vec(self.phi, Nvol)

        # Apply concentration and potential boundary conditions
        # Elyte value on the left and no-gradients on the right
        ctmp = np.hstack((self.portInLyte.c_lyte(), cvec, cvec[-1]))
        phitmp = np.hstack((self.portInLyte.phi_lyte(), phivec, phivec[-1]))

        Nm_edges, i_edges = get_interface_internal_fluxes(ctmp, phitmp, disc, config)

        disc["dxvec"][:] = 1.

        dvgNm = np.diff(Nm_edges) / disc["dxvec"]
        dvgi = np.diff(i_edges) / disc["dxvec"]

        for vInd in range(Nvol):
            # Mass Conservation (done with the anion, although "c" is neutral salt conc)
            eq = self.CreateEquation("interface_mass_cons_vol{vInd}".format(vInd=vInd))
            eq.Residual = disc["porosvec"][vInd]*dcdtvec[vInd] + (1./config["num"])*dvgNm[vInd]
            if config["interfaceModelType"] == "solid":
                eq.Residual += -config["kd"] * (config["cmax"] - cvec[vInd]) \
                    + config["kr"] * cvec[vInd] ** 2

            # Charge Conservation
            eq = self.CreateEquation("interface_charge_cons_vol{vInd}".format(vInd=vInd))
            eq.Residual = -dvgi[vInd]
            # Reaction out interface from last volume
            if vInd == Nvol - 1:
                # The volume of this particular particle
                Vj = config["psd_vol_FracVol"][self.trode][self.vInd,self.pInd]
                eq.Residual += config["zp"] * -(config["beta"][self.trode]
                                                * (1-config["poros"][self.trode])
                                                * config["P_L"][self.trode] * Vj
                                                * self.portInParticle.dcbardt())

            # Reaction entering the interface
            if vInd == 0:
                # The volume of this particular particle
                Vj = config["psd_vol_FracVol"][self.trode][self.vInd,self.pInd]
                eq.Residual -= config["zp"] * -(config["beta"][self.trode]
                                                * (1-config["poros"][self.trode])
                                                * config["P_L"][self.trode] * Vj
                                                * self.portInParticle.dcbardt())

        # last grid point of interface is output to particle
        eq = self.CreateEquation("c_interface_to_particle")
        eq.Residual = self.portOutInterfaceParticle.c_lyte() - ctmp[-1]

        eq = self.CreateEquation("phi_interface_to_particle")
        eq.Residual = self.portOutInterfaceParticle.phi_lyte() - phitmp[-1]

        eq = self.CreateEquation("Nm0_interface_to_elyte")
        eq.Residual = self.portOutInterfaceElyte.Nm0() - Nm_edges[0]

        eq = self.CreateEquation("i0_interface_to_elyte")
        eq.Residual = self.portOutInterfaceElyte.i0() - i_edges[1]

        for eq in self.Equations:
            eq.CheckUnitsConsistency = False


def get_interface_internal_fluxes(c, phi, disc, config):
    zp, zm, nup, num = config["zp"], config["zm"], config["nup"], config["num"]
    nu = nup + num
    T = config["T"]
    dxd1 = disc["dxd1"]
    eps_o_tau = disc["eps_o_tau"]

    # Get concentration at cell edges using weighted mean
    wt = utils.pad_vec(disc["dxvec"])
    c_edges_int = utils.weighted_linear_mean(c, wt)

    if config["interfaceModelType"] == "dilute":
        # Get porosity at cell edges using weighted harmonic mean
        eps_o_tau_edges = utils.weighted_linear_mean(eps_o_tau, wt)
        Dp = eps_o_tau_edges * config["Dp"]
        Dm = eps_o_tau_edges * config["Dm"]
#        Np_edges_int = nup*(-Dp*np.diff(c_lyte)/dxd1
#                            - Dp*zp*c_edges_int*np.diff(phi_lyte)/dxd1)
        Nm_edges_int = num*(-Dm*np.diff(c)/dxd1
                            - Dm/T*zm*c_edges_int*np.diff(phi)/dxd1)
        i_edges_int = (-((nup*zp*Dp + num*zm*Dm)*np.diff(c)/dxd1)
                       - (nup*zp**2*Dp + num*zm**2*Dm)/T*c_edges_int*np.diff(phi)/dxd1)
#        i_edges_int = zp*Np_edges_int + zm*Nm_edges_int
    elif config["interfaceModelType"] == "SM":
        D_fs, sigma_fs, thermFac, tp0 = getattr(props_elyte,config["interfaceSMset"])()[:-1]

        # Get diffusivity and conductivity at cell edges using weighted harmonic mean
        D_edges = utils.weighted_harmonic_mean(eps_o_tau*D_fs(c, T), wt)
        sigma_edges = utils.weighted_harmonic_mean(eps_o_tau*sigma_fs(c, T), wt)

        sp, n = config["sp"], config["n"]
        i_edges_int = -sigma_edges/T * (
            np.diff(phi)/dxd1
            + nu*T*(sp/(n*nup)+tp0(c_edges_int, T)/(zp*nup))
            * thermFac(c_edges_int, T)
            * np.diff(np.log(c))/dxd1
            )
        Nm_edges_int = num*(-D_edges*np.diff(c)/dxd1
                            + (1./(num*zm)*(1-tp0(c_edges_int, T))*i_edges_int))
    elif config["interfaceModelType"] == "solid":
        D_fs, sigma_fs, thermFac, tp0 = getattr(props_elyte, config["interfaceSMset"])()[:-1]

        # Get diffusivity at cell edges using weighted harmonic mean
        D_edges = utils.weighted_harmonic_mean(eps_o_tau * D_fs(c), wt)

        # sp, n = config["sp"], config["n_refTrode"]
        # D_fs is specified in solid_elyte_func in props_elyte.py
        Dm = config["Dm"]
        a_slyte = config["a_slyte"]
        k = 1.381e-23

        i_edges_int = (-((nup*zp*D_edges*(1/(1-c_edges_int)-a_slyte/(k*T)*2*c_edges_int)
                          + num*zm*Dm)*np.diff(c)/dxd1)
                       - (nup * zp ** 2 * D_edges + num * zm ** 2 * Dm) / T
                       * c_edges_int * np.diff(phi) / dxd1)
        Nm_edges_int = num * (-D_edges * np.diff(c) / dxd1
                              + (1. / (num * zm) * (1 - tp0(c_edges_int)) * i_edges_int))
    return Nm_edges_int, i_edges_int
