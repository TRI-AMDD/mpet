import numpy as np


def ICET(eta, c_sld, c_lyte, k0, E_A, T, act_R=None,
             act_lyte=None, lmbda=None, alpha=None):
    ecd = 0.03*(k0 * c_lyte**(1-alpha)
           * (1.0 - c_sld)*c_sld**alpha)
    Rate = ecd * np.exp(-E_A/T + E_A/1) * (np.exp(-alpha*eta/T) - np.exp((1-alpha)*eta/T))
    return Rate
