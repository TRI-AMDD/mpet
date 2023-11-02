import numpy as np
from .MHC_kfunc import MHC_kfunc


def ECIT_insul(eta, c_sld, c_lyte, k0, E_A, T, act_R=None,
         act_lyte=None, lmbda=None, alpha=None):
    # See Fraggedakis et al. 2020
    eta_f = eta + T*np.log(c_lyte/c_sld)
    ecd_extras = k0*(1-c_sld)
    Rrd = c_lyte*np.exp(-((lmbda+eta_f)**2)/(4*lmbda))
    Rox = c_sld*np.exp(-((lmbda-eta_f)**2)/(4*lmbda))
    Rate = np.exp(-E_A/T + E_A/1) * ecd_extras*(Rrd - Rox)
    return Rate
