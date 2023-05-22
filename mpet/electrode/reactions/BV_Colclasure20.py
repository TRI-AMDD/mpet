import numpy as np


def BV_Colclasure20(eta, c_sld, c_lyte, k0, E_A, T, act_R=None,
                    act_lyte=None, lmbda=None, alpha=None):
    """Implemented for the Finegan 2020/Colclasure 2020 model comparison
    for the NMC electrode"""
    ecd = k0 * (165*c_sld**5 - 752.4*c_sld**4 + 1241*c_sld**3
                - 941.7*c_sld**2 + 325*c_sld - 35.85) * (c_lyte/1.2)**0.5
    Rate = ecd * (np.exp(-alpha*eta/T) - np.exp((1-alpha)*eta/T))
    return Rate
