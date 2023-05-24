import numpy as np
from mpet.config import constants


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
        return (10**(D_00 + D_01 / (T - (T_g0 + T_g1 * Ce))
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
