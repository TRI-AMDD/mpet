import numpy as np
import daetools.pyDAE as dae
import scipy.special as spcl


def MHC_kfunc(eta, lmbda):
    a = 1. + np.sqrt(lmbda)
    if isinstance(eta, dae.pyCore.adouble):
        ERF = dae.Erf
    else:
        ERF = spcl.erf
    # evaluate with eta for oxidation, -eta for reduction
    return (np.sqrt(np.pi*lmbda) / (1 + np.exp(-eta))
            * (1. - ERF((lmbda - np.sqrt(a + eta**2))
               / (2*np.sqrt(lmbda)))))
