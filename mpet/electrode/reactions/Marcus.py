import numpy as np
import daetools.pyDAE as dae
from mpet.config.constants import reactions_epsilon as eps


def Marcus(eta, c_sld, c_lyte, k0, E_A, T, act_R=None,
           act_lyte=None, lmbda=None, alpha=None):
    if isinstance(c_sld, np.ndarray):
        c_sld = np.array([
            dae.Max(eps, c_sld[i]) for i in range(len(c_sld))])
    else:
        c_sld = dae.Max(eps, c_sld)
    alpha = 0.5*(1 + (T/lmbda) * np.log(dae.Max(eps, c_lyte)/c_sld))
    # We'll assume c_e = 1 (at the standard state for electrons)
#        ecd = ( k0 * np.exp(-lmbda/(4.*T)) *
#        ecd = ( k0 *
    ecd = (k0 * (1-c_sld)
           * c_lyte**((3-2*alpha)/4.)
           * c_sld**((1+2*alpha)/4.))
    Rate = (ecd * np.exp(-E_A/T + E_A/1) * np.exp(-eta**2/(4.*T*lmbda))
            * (np.exp(-alpha*eta/T) - np.exp((1-alpha)*eta/T)))
    return Rate
