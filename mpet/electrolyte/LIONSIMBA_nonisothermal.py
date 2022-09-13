from mpet.config import constants


def LIONSIMBA_nonisothermal():
    """ Set of parameters from LIONSIMBA validation. Torchio et al, 2016.
    """

    def tp0(c, T):
        return 0.364

    def sigma(c, T):
        c_dim = c*1000  # dimensionalized c
        T_dim = T*constants.T_ref
        r1 = -10.5
        r2 = 0.668e-3
        r3 = 0.494e-6
        r4 = 0.074
        r5 = -1.78e-5
        r6 = -8.86e-10
        r7 = -6.96e-5
        r8 = 2.8e-8
        sig_out = 1e-4 * c_dim * (r1 + r2*c_dim + r3*c_dim**2 + T_dim
                                  * (r4 + r5*c_dim + r6*c_dim**2)
                                  + T_dim**2 * (r7 + r8*c_dim))**2
        return sig_out  # m^2/s

    def D(c, T):
        c_dim = c*1000
        T_dim = T*constants.T_ref
        r1 = 4.43
        r2 = 54
        r3 = 229
        r4 = 5e-3
        r5 = 0.22e-3
        D_out = 1e-4 * 10**(-r1-r2/(T_dim-r3-r4*c_dim)-r5*c_dim)
        return D_out

    def therm_fac(c, T):
        return 1.

    Dref = D(constants.c_ref/1000, 1)

    def D_ndim(c, T):
        return D(c, T) / Dref

    def sigma_ndim(c, T):
        return sigma(c, T) * (
            constants.k*constants.T_ref/(constants.e**2*Dref*constants.N_A*(constants.c_ref)))
    return D_ndim, sigma_ndim, therm_fac, tp0, Dref
