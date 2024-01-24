import numpy as np


def BV_NMC_Park2021(eta, c_sld, c_lyte, k0, E_A, T, act_R=None,
                    act_lyte=None, lmbda=None, alpha=None):
    """Implemented for the Park, J 2021 model
    for the NMC electrode"""
    ecd = k0*act_lyte/(1+np.exp(25*(c_sld-0.6)))
    Rate = ecd * np.exp(-E_A/T + E_A/1) * (np.exp(-alpha*eta/T) - np.exp((1-alpha)*eta/T))
    return Rate
