import numpy as np

from mpet.config import constants


def Tesla_electrolyte():
    """ Set of parameters from Marc Berliner
    """

    def tp0(c, T):
        return 0.455

    def sigma(c, T):
        """ K_e evaluates the conductivity coefficients for the electrolyte phase [S/m]
        """
        # we use the dimensional form in the equations
        c = c*1000
        T = T*constants.T_ref
        return (-0.5182386306736273 +
                - 0.0065182740160006376 * c +
                + 0.0016958426698238335 * T +
                + 1.4464586693911396e-6 * c**2 +
                + 3.0336049598190174e-5 * c*T +
                + 3.046769609846814e-10 * c**3 +
                - 1.0493995729897995e-8 * c**2*T)

    def D(c, T):
        """
        D_e evaluates the diffusion coefficients for the electrolyte phase [m^2/s]
        """
        c = c*1000
        T = T*constants.T_ref
        a = np.array([1.8636977199162228e-8, -1.3917476882039536e-10, 3.1325506789441764e-14,
                      -7.300511963906146e-17, 5.119530992181613e-20, -1.1514201913496038e-23,
                      2.632651793626908e-13, -1.1262923552112963e-16, 2.614674626287283e-19,
                      -1.8321158900930459e-22, 4.110643579807774e-26])
        return a[0] + a[1]*T + a[2]*c*T + a[3]*c**2*T + a[4]*c**3*T + a[5]*c**4*T + \
            a[6]*T**2 + a[7]*c*T**2 + a[8]*c**2*T**2 + a[9]*c**3*T**2 + a[10]*c**4*T**2

    def therm_fac(c, T):
        return 1.

    Dref = D(constants.c_ref/1000, 1)

    def D_ndim(c, T):
        return D(c, T) / Dref

    def sigma_ndim(c, T):
        return sigma(c, T) * (
            constants.k*constants.T_ref/(constants.e**2*Dref*constants.N_A*(constants.c_ref)))
    return D_ndim, sigma_ndim, therm_fac, tp0, Dref
