import numpy as np


def BV_raw(eta, c_sld, c_lyte, k0, E_A, T, act_R=None,
           act_lyte=None, lmbda=None, alpha=None):
    ecd = k0
    Rate = ecd * np.exp(-E_A/T + E_A/1) * (np.exp(-alpha*eta/T) - np.exp((1-alpha)*eta/T))
    return Rate
