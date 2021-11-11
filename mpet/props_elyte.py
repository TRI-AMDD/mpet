r"""
This module provides functions defining properties of the ion-conducting
phase -- the electrolyte Manage functions for the parameters involved in
Stefan-Maxwell based concentrated electrolyte transport theory for
binary electrolytes.

Each electrolyte set must output functions for the following as a
function of c (electrolyte concentration, M)
 - Dchem [m^2/s] = the prefactor for grad(c) in species conservation
 - sigma [S/m] = the conductivity
 - (1 + dln(f_\pm)/dln(c)) = the "thermodynamic factor"
 - t_+^0 = the transference number of the cations
T in these equations is nondimensionalized wrt 298K
"""

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


def Doyle96_EC_DMC_1_2():
    """ Set of parameters from Doyle, Newman, et al. 1996.
    """
    def tp0(c, T):
        return 0.363

    def D(c, T):
        return 7.5e-11  # m^2/s

    def therm_fac(c, T):
        return 1.

    def sigma(c, T):
        r1 = 1.0793e-4
        r2 = 6.7461e-3
        r3 = 5.2245e-3
        r4 = 1.3605e-3
        r5 = 1.1724e-4
        k0 = r1 + r2*c - r3*c**2 + r4*c**3 - r5*c**4  # S/cm
        return(100*k0)

    Dref = D(constants.c_ref, 1)

    def D_ndim(c, T):
        return D(c, T) / Dref

    def sigma_ndim(c, T):
        return sigma(c, T) * (
            constants.k*constants.T_ref/(constants.e**2*Dref*constants.N_A*constants.c_ref))
    return D_ndim, sigma_ndim, therm_fac, tp0, Dref


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


def Colclasure20():
    def tp0(Ce, T):
        T = T*constants.T_ref  # T has units in this model
        Acoeff_Trfn = -0.0000002876102*T**2 + 0.0002077407*T - 0.03881203
        Bcoeff_Trfn = 0.000001161463*T**2 - 0.00086825*T + 0.1777266
        Ccoeff_Trfn = -0.0000006766258*T**2 + 0.0006389189*T + 0.3091761

        return (Acoeff_Trfn*Ce**2 + Bcoeff_Trfn*Ce
                + Ccoeff_Trfn)

    def D(Ce, T):
        T = T*constants.T_ref  # T has units in this model
        D_00 = -0.5688226
        D_01 = -1607.003
        D_10 = -0.8108721
        D_11 = 475.291
        D_20 = -0.005192312
        D_21 = -33.43827
        T_g0 = -24.83763
        T_g1 = 64.07366
        return(10**(D_00 + D_01 / (T - (T_g0 + T_g1 * Ce))
                    + (D_10 + D_11 / (T - (T_g0 + T_g1 * Ce)))*Ce
                    + (D_20 + D_21 / (T - (T_g0 + T_g1 * Ce)))*Ce**2)*0.0001)

    def therm_fac(c, T):
        T = T*constants.T_ref  # T has units in this model
        return 1 + 0.341*np.exp(261/T)-0.00225*np.exp(1360/T)*c + 0.54*np.exp(329/T)*c**2

    def sigma(Ce, T):
        T = T*constants.T_ref  # T has units in this model
        Acoeff_Kappa = 0.0001909446*T**2 - 0.08038545*T + 9.00341
        Bcoeff_Kappa = -0.00000002887587*T**4 + 0.00003483638*T**3 - 0.01583677*T**2 \
            + 3.195295*T - 241.4638
        Ccoeff_Kappa = 0.00000001653786*T**4 - 0.0000199876*T**3 + 0.009071155*T**2 \
            - 1.828064*T + 138.0976
        Dcoeff_Kappa = -0.000000002791965*T**4 + 0.000003377143 * \
            T**3 - 0.001532707*T**2 + 0.3090003*T - 23.35671

        return (Ce * (Acoeff_Kappa + Bcoeff_Kappa * Ce
                      + Ccoeff_Kappa * Ce**2 + Dcoeff_Kappa * Ce**3))

    Dref = D(constants.c_ref/1000, 1)

    def D_ndim(c, T):
        return D(c, T) / Dref

    def sigma_ndim(c, T):
        return sigma(c, T) * (
            constants.k*constants.T_ref/(constants.e**2*Dref*constants.N_A*constants.c_ref))
    return D_ndim, sigma_ndim, therm_fac, tp0, Dref


def solid_elyte():
    """
    Solid Electrolyte version, several sources for different params
    """
    # LCO: kappa_ndim is not actually returned, so unclear
    # what it should be used for
    # related functions and values commented out for now

    tp0 = 0.9  # tp0 is constant but a callable must be returned
    D = 1.19e-11  # m^2/s
    # kappa = 1.2e-3  # S/m

    # def kappa_valoen_reimers(c):
    #     k00, k01, k02, k10, k11, k12, k20, k21 = (-10.5, 0.0740, -6.96e-5, 0.668,
    #                                               -0.0178, 2.80e-5, 0.494, -8.86e-4)
    #     out = c * (k00 + k01*Tref + k02*Tref**2
    #                + k10*c + k11*c*Tref + k12*c*Tref**2
    #                + k20*c**2 + k21*c**2*Tref)**2  # mS/cm
    #     out *= 0.1  # S/m
    #     return out

    Ign1, sigma_ndim, thermFac, Ign3, Dref = valoen_reimers()
    D_ndim = D / Dref
    # kappa_ndim = lambda c: kappa / kappa_valoen_reimers(c)

    # D_ndim and tp0 are constants, but a callable must be returned
    return lambda c: D_ndim, sigma_ndim, thermFac, lambda c: tp0, Dref


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
