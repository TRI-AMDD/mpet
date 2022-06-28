from mpet.config import constants


def LIONSIMBA_isothermal():
    """ Set of parameters from LIONSIMBA validation. Torchio et al, 2016.
    """
    def tp0(c, T):
        return 0.364

    def sigma(c, T):
        ce = c*1000  # dimensionalized c
        return (4.1253e-2 + 5.007e-4*ce - 4.7212e-7*ce**2
                + 1.5094e-10*ce**3 - 1.6018*1e-14*ce**4)  # S/m

    def D(c, T):
        return 7.5e-10  # m^2/s
    # isothermal at 298 K

    def therm_fac(c, T):
        return 1.

    Dref = D(constants.c_ref/1000, 1)

    def D_ndim(c, T):
        return D(c, T) / Dref

    def sigma_ndim(c, T):
        return sigma(c, T) * (
            constants.k*constants.T_ref/(constants.e**2*Dref*constants.N_A*(constants.c_ref)))
    return D_ndim, sigma_ndim, therm_fac, tp0, Dref
