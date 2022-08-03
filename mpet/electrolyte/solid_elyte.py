from .valoen_reimers import valoen_reimers


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
