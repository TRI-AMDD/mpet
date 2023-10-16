import numpy as np


def BV_NMC811(eta, c_sld, c_lyte, k0, E_A, T, act_R=None,
                    act_lyte=None, lmbda=None, alpha=None):
    """NMC811 Butler-Volmer reaction rate """
    # goes from 1 to 6
    # has to be multiplied by k0 in A/m^2
    # ecd = k0*(16.69*c_sld**5 + -279.99*c_sld**4
    #            + 687.19*c_sld**3 + -634.2*c_sld**2
    #              + 234.07*c_sld + -23.92) * (c_lyte)**0.5
    ecd = k0*(87.15*c_sld**5 + -445.05*c_sld**4 
              + 815.43*c_sld**3 + -669.95*c_sld**2 
              + 236.09*c_sld + -23.67) * (c_lyte)**0.5
    Rate = ecd * (np.exp(-alpha*eta/T) - np.exp((1-alpha)*eta/T))
    return Rate
