from mpet.config import constants


def valoen_reimers():
    """ Set of parameters from Valoen and Reimers 2005 """
    def tp0(c, T):
        return 0.38

    def D(c, T):
        return (
            10**(-4) * 10**(-4.43 - 54/((T*constants.T_ref) - (229 + 5*c)) - 0.22*c))  # m^2/s

    def therm_fac(c, T):
        tmp = 0.601 - 0.24*c**(0.5) + 0.982*(1 - 0.0052*((T*constants.T_ref) - 294))*c**(1.5)
        return tmp/(1-tp0(c, T))

    def sigma(c, T):
        (k00, k01, k02,
         k10, k11, k12,
         k20, k21) = (
             -10.5, 0.0740, -6.96e-5,
             0.668, -0.0178, 2.80e-5,
             0.494, -8.86e-4)
        out = c * (k00 + k01*(T*constants.T_ref) + k02*(T*constants.T_ref)**2
                   + k10*c + k11*c*(T*constants.T_ref) + k12*c*(T*constants.T_ref)**2
                   + k20*c**2 + k21*c**2*(T*constants.T_ref))**2  # mS/cm
        out *= 0.1  # S/m
        return out
    Dref = D(constants.c_ref/1000, 1)

    def D_ndim(c, T):
        return D(c, T) / Dref

    def sigma_ndim(c, T):
        return sigma(c, T) * (
            constants.k*constants.T_ref/(constants.e**2*Dref*constants.N_A*(constants.c_ref)))
    return D_ndim, sigma_ndim, therm_fac, tp0, Dref
