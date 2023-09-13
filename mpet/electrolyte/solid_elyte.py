from .valoen_reimers import valoen_reimers


def solid_elyte():
    """
    Solid Electrolyte version
    """
    tp0 = 0.99999
    D = 1.19e-11  # m^2/s

    Ign1, sigma_ndim, thermFac, Ign3, Dref = valoen_reimers()
    D_ndim = D / Dref

    # D_ndim and tp0 are constants, but a callable must be returned
    return lambda c: D_ndim, sigma_ndim, thermFac, tp0, Dref
