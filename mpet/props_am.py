"""This module handles properties associated with the active materials.
Only helper functions are defined here.
Diffusion functions are defined in mpet.electrode.diffusion
Chemical potential functions are defined in mpet.electrode.materials"""
import types
import numpy as np

import mpet.geometry as geo
from mpet.config import constants
from mpet.utils import import_function


class muRfuncs():
    """ This class defines functions which describe the chemical
    potential of active materials.
    Each function describes a particular material. The
    chemical potential is directly related to the open circuit voltage
    (OCV) for solid solution materials.
    In each function, muR_ref is the offset (non-dimensional) chemical
    potential by which the returned value will be shifted (for
    numerical convenience).
    Each material function returns:
        muR -- chemical potential
        actR -- activity (if applicable, else None)
    """

    def __init__(self, config, trode, ind=None):
        """config is the full dictionary of
        parameters for the electrode particles, as made for the
        simulations. trode is the selected electrode.
        ind is optinally the selected particle, provided as (vInd, pInd)
        """
        self.config = config
        self.trode = trode
        self.ind = ind
        # eokT and kToe are the reference values for scalings
        self.eokT = constants.e / (constants.k * constants.T_ref)
        self.kToe = 1. / self.eokT

        # If the user provided a filename with muRfuncs, try to load
        # the function from there, otherwise load it from this class
        filename = self.get_trode_param("muRfunc_filename")
        muRfunc_name = self.get_trode_param("muRfunc")
        if filename is None:
            # the function will be loaded from the materials folder
            muRfunc = import_function(None, muRfunc_name,
                                      mpet_module=f"mpet.electrode.materials.{muRfunc_name}")
        else:
            muRfunc = import_function(filename, muRfunc_name)

        # We have to make sure the function knows what 'self' is with
        # the types.MethodType function
        self.muRfunc = types.MethodType(muRfunc, self)

    def get_trode_param(self, item):
        """
        Shorthand to retrieve electrode-specific value
        """
        value = self.config[self.trode, item]
        # check if it is a particle-specific parameter
        if self.ind is not None and item in self.config.params_per_particle:
            value = value[self.ind]
        return value

    def get_muR_from_OCV(self, OCV, muR_ref):
        return -self.eokT*OCV + muR_ref

    ######
    # Helper functions
    ######

    def ideal_sln(self, y, T):
        """ Helper function: Should not be called directly from
        simulation. Call a specific material instead. """
        muR = T*np.log(y/(1-y))
        return muR

    def reg_sln(self, y, T, Omga):
        """ Helper function """
        muR_IS = self.ideal_sln(y, T)
        enthalpyTerm = Omga*(1-2*y)
        muR = muR_IS + enthalpyTerm
        return muR

    def graphite_2param_homog(self, y, T, Omga, Omgb, Omgc, EvdW):
        """ Helper function """
        y1, y2 = y
        muR1 = self.reg_sln(y1, T, Omga)
        muR2 = self.reg_sln(y2, T, Omga)
        muR1 += Omgb*y2 + Omgc*y2*(1-y2)*(1-2*y1)
        muR2 += Omgb*y1 + Omgc*y1*(1-y1)*(1-2*y2)
        muR1 += EvdW * (30 * y1**2 * (1-y1)**2)
        muR2 += EvdW * (30 * y2**2 * (1-y2)**2)
        return (muR1, muR2)

    def graphite_1param_homog(self, y, T, Omga, Omgb):
        """ Helper function """
        width = 5e-2
        tailScl = 5e-2
        muLtail = -tailScl*1./(y**(0.85))
        muRtail = tailScl*1./((1-y)**(0.85))
        slpScl = 0.45
        muLlin = slpScl*Omga*4*(0.26-y)*step_down(y, 0.5, width)
        muRlin = (slpScl*Omga*4*(0.74-y) + Omgb)*step_up(y, 0.5, width)
        muR = muLtail + muRtail + muLlin + muRlin
        return muR

    def graphite_1param_homog_2(self, y, T, Omga, Omgb):
        """ Helper function """
        width = 5e-2
        tailScl = 5e-2
        slpScl = 0.45
        muLtail = -tailScl*1./(y**(0.85))
        muRtail = tailScl*1./((1-y)**(0.85))
        muLlin = (slpScl*Omga*12*(0.40-y)
                  * step_down(y, 0.49, 0.9*width)*step_up(y, 0.35, width))
        muRlin = (slpScl*Omga*4*(0.74-y) + Omgb)*step_up(y, 0.5, width)
        muLMod = (0.
                  + 40*(-np.exp(-y/0.015))
                  + 0.75*(np.tanh((y-0.17)/0.02) - 1)
                  + 1.0*(np.tanh((y-0.22)/0.040) - 1)
                  )*step_down(y, 0.35, width)
        muR = muLMod + muLtail + muRtail + muLlin + muRlin
        return muR

    def graphite_1param_homog_3(self, y, T, Omga, Omgb):
        """ Helper function with low hysteresis and soft tail """
        width = 5e-2
        tailScl = 5e-2
        muLtail = -tailScl*1./(y**(0.85))
        muRtail = tailScl*1./((1-y)**(0.85))
        muRtail = 1.0e1*step_up(y, 1.0, 0.045)
        muLlin = (0.15*Omga*12*(0.40-y**0.98)
                  * step_down(y, 0.49, 0.9*width)*step_up(y, 0.35, width))
        muRlin = (0.1*Omga*4*(0.74-y) + 0.90*Omgb)*step_up(y, 0.5, 0.4*width)
        muLMod = (0.
                  + 40*(-np.exp(-y/0.015))
                  + 0.75*(np.tanh((y-0.17)/0.02) - 1)
                  + 1.0*(np.tanh((y-0.22)/0.040) - 1)
                  )*step_down(y, 0.35, width)
        muR = 0.18 + muLMod + muLtail + muRtail + muLlin + muRlin
        return muR
    
    def reg_sol_red(self, y, T, Omga, theta):
        tail = 0.075
        delta_y = 0.1
        m = (4 - 2*Omga)
        muR = (
            (theta - (m/2)+ m*y)*step_down(y, (1+delta_y), tail)
            + theta*step_down(y, (0-delta_y), tail)
        )
        return muR

    def LMFP_red(self, y, T, Omga, Omgb, theta_1, theta_2, stoich):

        Omegc = (Omga + Omgb)*0.5
        m_1 = (4 - 2*Omga*stoich)
        m_2 = (4 - 2*Omgb*(1-stoich))
        tail_out = 0.015
        tail_mid = 0.015

        # theta_1 = -4.09
        # theta_2 = -3.422
        mu_reduced = (
            (theta_1-(m_1*stoich/2)+m_1*y)*step_down(y,stoich,tail_mid)
            + (theta_2-(m_2*(1-stoich)/2)-m_2*stoich+m_2*y)*step_up(y,stoich,tail_mid)
            + theta_1*step_down(y,-0.01,tail_out)
            + theta_2*step_down(y,1.01,tail_out) - theta_2
            - stoich*Omegc*step_up(y,stoich,tail_mid)
            + (1-stoich)*Omegc*step_down(y,stoich,tail_mid)
            )
        return mu_reduced

    def non_homog_rect_fixed_csurf(self, y, ybar, B, kappa, ywet):
        """ Helper function """
        if isinstance(y, np.ndarray):
            N = len(y)
            ytmp = np.empty(N+2, dtype=object)
            ytmp[1:-1] = y
            ytmp[0] = ywet
            ytmp[-1] = ywet
            dxs = 1./N
            curv = np.diff(ytmp, 2)/(dxs**2)
            muR_nh = -kappa*curv + B*(y - ybar)
        elif (isinstance(y, tuple) and len(y) == 2
                and isinstance(y[0], np.ndarray)):
            stoich_1 = self.get_trode_param("stoich_1")
            stoich_2 = 1 - stoich_1
            ybar_avg = stoich_1*ybar[0]+stoich_2*ybar[1]
            y_avg = stoich_1*y[0]+stoich_2*y[1]
            N = len(y[0])
            kappa1 = kappa[0]
            kappa2 = kappa[1]
            B1 = B[0]
            B2 = B[1]
            B_avg = stoich_1*B1 + stoich_2*B2
            ytmp1 = np.empty(N+2, dtype=object)
            ytmp1[1:-1] = y[0]
            ytmp1[0] = ywet
            ytmp1[-1] = ywet
            dxs = 1./N
            curv1 = np.diff(ytmp1, 2)/(dxs**2)
            muR1_nh = -stoich_1*kappa1*curv1 + B_avg*(y_avg - ybar_avg)
            ytmp2 = np.empty(N+2, dtype=object)
            ytmp2[1:-1] = y[1]
            ytmp2[0] = ywet
            ytmp2[-1] = ywet
            curv2 = np.diff(ytmp2, 2)/(dxs**2)
            muR2_nh = -stoich_2*kappa2*curv2 + B_avg*(y_avg - ybar_avg)

            muR_nh = (muR1_nh, muR2_nh)
        return muR_nh

    # def non_homog_round_wetting(self, y, ybar, B, kappa, beta_s, shape, r_vec):
    #         """ Helper function """
    #         dr = r_vec[1] - r_vec[0]
    #         Rs = 1.
    #         curv = geo.calc_curv(y, dr, r_vec, Rs, beta_s, shape)
    #         muR_nh = B*(y - ybar) - kappa*curv
    #         return muR_nh

    def non_homog_round_wetting(self, y, ybar, B, kappa, beta_s, shape, r_vec):
        """ Helper function """
        mod1var, mod2var = False, False
        if isinstance(y, np.ndarray):
            mod1var = True
        elif (isinstance(y, tuple) and len(y) == 2
                and isinstance(y[0], np.ndarray)):
            mod2var = True
        if mod1var:
            dr = r_vec[1] - r_vec[0]
            Rs = 1.
            curv = geo.calc_curv(y, dr, r_vec, Rs, beta_s, shape)
            muR_nh = B*(y - ybar) - kappa*curv
        elif mod2var:
            if self.get_trode_param("stoich_1") is not None:
                stoich_1 = self.get_trode_param("stoich_1")
                stoich_2 = 1 - stoich_1
            else:
                stoich_1 = stoich_2 = 0.5
            y1 = y[0]
            y2 = y[1]
            kappa1 = kappa[0]
            kappa2 = kappa[1]
            dr = r_vec[1] - r_vec[0]
            Rs = 1.
            # wetting needs to be fixed
            curv1 = geo.calc_curv(y1, dr, r_vec, Rs, beta_s, shape)
            curv2 = geo.calc_curv(y2, dr, r_vec, Rs, beta_s, shape)
            B1 = B[0]
            B2 = B[1]
            y_avg = stoich_1*y1+stoich_2*y2
            ybar_avg = stoich_1*ybar[0]+stoich_2*ybar[1]
            B_avg = stoich_1*B1 + stoich_2*B2

            muR1_nh = B_avg*(y_avg - ybar_avg) - stoich_1*kappa1*curv1
            muR2_nh = B_avg*(y_avg - ybar_avg) - stoich_2*kappa2*curv2

            muR_nh = (muR1_nh, muR2_nh)
        return muR_nh

    def general_non_homog(self, y, ybar):
        """ Helper function """
        ptype = self.get_trode_param("type")
        mod1var, mod2var = False, False
        if isinstance(y, np.ndarray):
            mod1var = True
            N = len(y)
        elif (isinstance(y, tuple) and len(y) == 2
                and isinstance(y[0], np.ndarray)):
            mod2var = True
            N = len(y[0])
        else:
            raise Exception("Unknown input type")
        if ("homog" not in ptype) and (N > 1):
            shape = self.get_trode_param("shape")
            kappa = self.get_trode_param("kappa")
            B = self.get_trode_param("B")
            if shape == "C3":
                if mod1var:
                    cwet = self.get_trode_param("cwet")
                    muR_nh = self.non_homog_rect_fixed_csurf(
                        y, ybar, B, kappa, cwet)
                elif mod2var:
                    # only type using this is ACR2
                    if self.get_trode_param("kappa1") is not None:
                        kappa1 = self.get_trode_param("kappa1")
                        kappa2 = self.get_trode_param("kappa2")
                    else:
                        kappa1 = kappa2 = self.get_trode_param("kappa")
                    kappa = (kappa1,kappa2)
                    if self.get_trode_param("B1") is not None:
                        B1 = self.get_trode_param("B1")
                        B2 = self.get_trode_param("B2")
                    else:
                        B1 = B2 = self.get_trode_param("B")
                    B = (B1,B2)
                    cwet = self.get_trode_param("cwet")
                    muR_nh = self.non_homog_rect_fixed_csurf(
                        y, ybar, B, kappa, cwet)
            # elif shape in ["cylinder", "sphere"]:
            #     beta_s = self.get_trode_param("beta_s")
            #     r_vec = geo.get_unit_solid_discr(shape, N)[0]
            #     if mod1var:
            #         muR_nh = self.non_homog_round_wetting(
            #             y, ybar, B, kappa, beta_s, shape, r_vec)
            #     elif mod2var:
            #         muR1_nh = self.non_homog_round_wetting(
            #             y[0], ybar[0], B, kappa, beta_s, shape, r_vec)
            #         muR2_nh = self.non_homog_round_wetting(
            #             y[1], ybar[1], B, kappa, beta_s, shape, r_vec)
            #         muR_nh = (muR1_nh, muR2_nh)
            elif shape in ["cylinder", "sphere"]:
                beta_s = self.get_trode_param("beta_s")
                r_vec = geo.get_unit_solid_discr(shape, N)[0]
                if mod1var:
                    muR_nh = self.non_homog_round_wetting(
                        y, ybar, B, kappa, beta_s, shape, r_vec)
                elif mod2var:
                    if self.get_trode_param("kappa1") is not None:
                        kappa1 = self.get_trode_param("kappa1")
                        kappa2 = self.get_trode_param("kappa2")
                    else:
                        kappa1 = kappa2 = self.get_trode_param("kappa")
                    kappa = (kappa1,kappa2)
                    if self.get_trode_param("B1") is not None:
                        B1 = self.get_trode_param("B1")
                        B2 = self.get_trode_param("B2")
                    else:
                        B1 = B2 = self.get_trode_param("B")
                    B = (B1,B2)
                    muR_nh = self.non_homog_round_wetting(
                        y, ybar, B, kappa, beta_s, shape, r_vec)
        else:  # homogeneous particle
            if mod1var:
                muR_nh = 0*y
            elif mod2var:
                muR_nh = (0*y[0], 0*y[1])
        return muR_nh


def step_down(x, xc, delta):
    return 0.5*(-np.tanh((x - xc)/delta) + 1)


def step_up(x, xc, delta):
    return 0.5*(np.tanh((x - xc)/delta) + 1)
