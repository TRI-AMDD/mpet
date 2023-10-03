import numpy as np

import daetools.pyDAE as dae
from mpet import ports, utils
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
        T = self.config["T"]

        # Concentration: if elyte solid not continuous over
        # interface if liquid continous over interface.
        if config["interfaceModelType"] == "solid" or config["elyteModelType"] == "solid":
            ctmp = np.hstack((cvec[0], cvec, cvec[-1]))
        else:
            ctmp = np.hstack((self.portInLyte.c_lyte(), cvec, cvec[-1]))

        # Electrical potential
        if config["elyteModelType"] == "dilute":
            phitmp = np.hstack((self.portInLyte.phi_lyte()
                                + T * np.log(self.portInLyte.c_lyte()), phivec, phivec[-1]))
        elif config["interfaceModelType"] == "dilute":
            phitmp = np.hstack((self.portInLyte.phi_lyte()
                                - T * np.log(cvec[0]), phivec, phivec[-1]))
        else:
            phitmp = np.hstack((self.portInLyte.phi_lyte(), phivec, phivec[-1]))

        Nm_edges, i_edges = get_interface_internal_fluxes(ctmp, phitmp, disc, config)

        # The reaction rate per volume (Rvp) is normalized to the total length of the electrode.
        dlc = config["L"][self.trode]/config["Nvol"][self.trode]
        disc["dxvec"][:] = 1

        dvgNm = np.diff(Nm_edges) / disc["dxvec"]
        dvgi = np.diff(i_edges) / disc["dxvec"]

        for vInd in range(Nvol):
            # Mass Conservation (done with the anion, although "c" is neutral salt conc)
            eq = self.CreateEquation("interface_mass_cons_vol{vInd}".format(vInd=vInd))
            eq.Residual = disc["porosvec"][vInd]*dcdtvec[vInd] + (1./config["num"])*dvgNm[vInd]

            # Charge Conservation
            eq = self.CreateEquation("interface_charge_cons_vol{vInd}".format(vInd=vInd))
            eq.Residual = -dvgi[vInd]
            # Reaction out interface from last volume
            if vInd == Nvol - 1:
                # The volume of this particular particle
                Vj = config["psd_vol_FracVol"][self.trode][self.vInd,self.pInd]
                eq.Residual += dlc * config["zp"] * -(config["beta"][self.trode]
                                                      * (1-config["poros"][self.trode])
                                                      * config["P_L"][self.trode] * Vj
                                                      * self.portInParticle.dcbardt())

            # Reaction entering the interface
            if vInd == 0:
                # The volume of this particular particle
                Vj = config["psd_vol_FracVol"][self.trode][self.vInd,self.pInd]
                eq.Residual -= dlc * config["zp"] * -(config["beta"][self.trode]
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
        # eq.Residual = self.portOutInterfaceElyte.i0() - i_edges[1]
        eq.Residual = self.portOutInterfaceElyte.i0() - i_edges[1] / dlc

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
        Dp = eps_o_tau_edges * config["Dp_i"]
        Dm = eps_o_tau_edges * config["Dm_i"]
#        Np_edges_int = nup*(-Dp*np.diff(c_lyte)/dxd1
#                            - Dp*zp*c_edges_int*np.diff(phi_lyte)/dxd1)
        Nm_edges_int = num*(-Dm*np.diff(c)/dxd1
                            - Dm/T*zm*c_edges_int*np.diff(phi)/dxd1)
        i_edges_int = (-((nup*zp*Dp + num*zm*Dm)*np.diff(c)/dxd1)
                       - (nup*zp**2*Dp + num*zm**2*Dm)/T*c_edges_int*np.diff(phi)/dxd1)
#        i_edges_int = zp*Np_edges_int + zm*Nm_edges_int
    elif config["interfaceModelType"] == "SM":
        SMset = config["SMset"]
        elyte_function = utils.import_function(config["SMset_filename"], SMset,
                                               mpet_module=f"mpet.electrolyte.{SMset}")
        D_fs, sigma_fs, thermFac, tp0 = elyte_function()[:-1]

        # Get diffusivity and conductivity at cell edges using weighted harmonic mean
        D_edges = utils.weighted_harmonic_mean(eps_o_tau*D_fs(c, T), wt)
        sigma_edges = utils.weighted_harmonic_mean(eps_o_tau*sigma_fs(c, T), wt)

        sp, n = config["sp"], config["n"]
        # there is an error in the MPET paper, temperature dependence should be
        # in sigma and not outside of sigma
        i_edges_int = -sigma_edges * (
            np.diff(phi)/dxd1
            + nu*T*(sp/(n*nup)+tp0(c_edges_int, T)/(zp*nup))
            * thermFac(c_edges_int, T)
            * np.diff(np.log(c))/dxd1
            )
        Nm_edges_int = num*(-D_edges*np.diff(c)/dxd1
                            + (1./(num*zm)*(1-tp0(c_edges_int, T))*i_edges_int))

    elif config["interfaceModelType"] == "solid":
        SMset = config["SMset"]
        elyte_function = utils.import_function(config["SMset_filename"], SMset,
                                               mpet_module=f"mpet.electrolyte.{SMset}")
        D_fs, sigma_fs, thermFac, tp0 = elyte_function()[:-1]

        a_slyte = config["a_slyte"]
        tp0 = 0.99999

        c_edges_int_norm = c_edges_int / config["cmax_i"]

        # Get diffusivity at cell edges using weighted harmonic mean
        # D_edges = utils.weighted_harmonic_mean(eps_o_tau * D_fs(c_lyte), wt)
        eps_o_tau_edges = utils.weighted_linear_mean(eps_o_tau, wt)
        # sp, n = ndD["sp"], ndD["n_refTrode"]
        # D_fs is specified in solid_elyte_func in props_elyte.py
        Dp = eps_o_tau_edges * config["Dp_i"]
        Dm = (zp * Dp - zp * Dp * tp0) / (tp0 * zm)

        Dp0 = Dp / (1-c_edges_int_norm)  # should be c0/cmax

        Dchemp = Dp0 * (1 - 2 * a_slyte * c_edges_int_norm + 2 * a_slyte * c_edges_int_norm**2)
        Dchemm = Dm

        Damb = (zp * Dp * Dchemm + zm * Dm * Dchemp) / (zp * Dp - zm * Dm)

        i_edges_int = (-((nup*zp*Dchemp + num*zm*Dchemm)*np.diff(c)/dxd1)
                       - (nup * zp ** 2 * Dp0 * (1 - c_edges_int_norm) + num * zm ** 2 * Dm) / T
                       * c_edges_int * np.diff(phi) / dxd1)

        Nm_edges_int = num * (-Damb * np.diff(c) / dxd1
                              + (1. / (num * zm) * (1 - tp0) * i_edges_int))
    return Nm_edges_int, i_edges_int
