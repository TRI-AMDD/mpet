from mpet.config import constants


def Doyle96_EC_DMC_2_1():
    """ Set of parameters from Doyle, Newman, et al. 1996.
    """
    def tp0(c, T):
        return 0.363

    def D(c, T):
        return 7.5e-11  # m^2/s

    def therm_fac(c, T):
        return 1.

    def sigma(c, T):
        r1 = 4.1253e-4
        r2 = 5.007e-3
        r3 = 4.7212e-3
        r4 = 1.5094e-3
        r5 = 1.6018e-4
        k0 = r1 + r2*c - r3*c**2 + r4*c**3 - r5*c**4  # S/cm
        return(100*k0)

    Dref = D(constants.c_ref, 1)

    def D_ndim(c, T):
        return D(c, T) / Dref

    def sigma_ndim(c, T):
        return sigma(c, T) * (
            constants.k*constants.T_ref/(constants.e**2*Dref*constants.N_A*constants.c_ref))
    return D_ndim, sigma_ndim, therm_fac, tp0, Dref
