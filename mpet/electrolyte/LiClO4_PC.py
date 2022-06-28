import numpy as np
from mpet.config import constants


def LiClO4_PC():
    """ Set of parameters from Fuller, Doyle, Newman 1994, with
    conductivity directly from dualfoil5.2.f
    """
    def tp0(c, T):
        return 0.2

    def D(c, T):
        return 2.58e-10  # m^2/s

    def therm_fac(c, T):
        # therm_fac adjusted to account for missing a factor of 2 in Eq A-2
        return 0.5

    def sigma(cin, T):
        c = cin * 1000  # mol/m^3
        p_max = 0.542
        p_u = 0.6616
        a = 0.855
        b = -0.08
        rho = 1.2041e3
        out = 0.0001 + c**a * (
            p_max*(1./(rho*p_u))**a
            * np.exp(b*(c/rho - p_u)**2
                     - (a/p_u)*(c/rho - p_u)))  # S/m
        return out
    Dref = D(constants.c_ref/1000, 1)

    def D_ndim(c, T):
        return D(c, T) / Dref

    def sigma_ndim(c, T):
        return sigma(c, T) * (
            constants.k*constants.T_ref/(constants.e**2*Dref*constants.N_A*(constants.c_ref)))
    return D_ndim, sigma_ndim, therm_fac, tp0, Dref
