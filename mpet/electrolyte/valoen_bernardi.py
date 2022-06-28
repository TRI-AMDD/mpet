from mpet.config import constants
from .valoen_reimers import valoen_reimers


def valoen_bernardi():
    """ Set of parameters from Bernardi and Go 2011, indirectly from
    Valoen and Reimers 2005. The only change from Valoen and Reimers
    is the conductivity.
    """
    D_ndim, Ign, therm_fac, tp0, Dref = valoen_reimers()

    def sigma(c, T):
        (k00, k01, k02,
         k10, k11, k12,
         k20, k21) = (
            -8.2488, 0.053248, -0.000029871,
            0.26235, -0.0093063, 0.000008069,
            0.22002, -0.0001765)
        out = c * (k00 + k01*(T*constants.T_ref) + k02*(T*constants.T_ref)**2
                   + k10*c + k11*c*(T*constants.T_ref) + k12*c*(T*constants.T_ref)**2
                   + k20*c**2 + k21*c**2*(T*constants.T_ref))**2  # mS/cm
        out *= 0.1  # S/m
        return out

    def sigma_ndim(c, T):
        return sigma(c, T) * (
            constants.k*constants.T_ref/(constants.e**2*Dref*constants.N_A*(constants.c_ref)))
    return D_ndim, sigma_ndim, therm_fac, tp0, Dref
