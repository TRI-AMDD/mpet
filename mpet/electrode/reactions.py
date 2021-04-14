import daetools.pyDAE as dae
import numpy as np
import scipy.special as spcl

eps = -1e-12

# Reaction rate functions
# Should have the same name as the rxnType config field


def BV(eta, c_sld, c_lyte, k0, E_A, T, act_R=None,
       act_lyte=None, lmbda=None, alpha=None):
    if act_R is None:
        act_R = c_sld/(1-c_sld)
    gamma_ts = (1./(1-c_sld))
    ecd = (k0 * act_lyte**(1-alpha)
           * act_R**(alpha) / gamma_ts)
    Rate = ecd * np.exp(-E_A/T + E_A/1) * (np.exp(-alpha*eta/T) - np.exp((1-alpha)*eta/T))
    return Rate


def BV_gMod01(eta, c_sld, c_lyte, k0, E_A, T, act_R=None,
              act_lyte=None, lmbda=None, alpha=None):
    if act_R is None:
        act_R = c_sld/(1-c_sld)
    gamma_ts = (1./(c_sld*(1-c_sld)))
    ecd = (k0 * act_lyte**(1-alpha)
           * act_R**(alpha) / gamma_ts)
    Rate = ecd * np.exp(-E_A/T + E_A/1) * (np.exp(-alpha*eta/T) - np.exp((1-alpha)*eta/T))
    return Rate


def BV_raw(eta, c_sld, c_lyte, k0, E_A, T, act_R=None,
           act_lyte=None, lmbda=None, alpha=None):
    ecd = k0
    Rate = ecd * np.exp(-E_A/T + E_A/1) * (np.exp(-alpha*eta/T) - np.exp((1-alpha)*eta/T))
    return Rate


def BV_mod01(eta, c_sld, c_lyte, k0, E_A, T, act_R=None,
             act_lyte=None, lmbda=None, alpha=None):
    ecd = (k0 * c_lyte**(1-alpha)
           * (1.0 - c_sld)**(1 - alpha) * c_sld**alpha)
    Rate = ecd * np.exp(-E_A/T + E_A/1) * (np.exp(-alpha*eta/T) - np.exp((1-alpha)*eta/T))
    return Rate


def BV_mod02(eta, c_sld, c_lyte, k0, E_A, T, act_R=None,
             act_lyte=None, lmbda=None, alpha=None):
    ecd = (k0 * c_lyte**(1-alpha)
           * (0.5 - c_sld)**(1 - alpha) * c_sld**alpha)
    Rate = ecd * np.exp(-E_A/T + E_A/1) * (np.exp(-alpha*eta/T) - np.exp((1-alpha)*eta/T))
    return Rate


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


def MHC(eta, c_sld, c_lyte, k0, E_A, T, act_R=None,
        act_lyte=None, lmbda=None, alpha=None):
    # See Zeng, Smith, Bai, Bazant 2014
    # Convert to "MHC overpotential"
    k0 = k0/MHC_kfunc(0., lmbda)
    eta_f = eta + T*np.log(c_lyte/c_sld)
    gamma_ts = 1./(1. - c_sld)
    alpha = 0.5
    ecd_extras = act_lyte**(1-alpha) * act_R**(alpha) / (gamma_ts*np.sqrt(c_lyte*c_sld))
    if isinstance(eta, np.ndarray):
        Rate = np.empty(len(eta), dtype=object)
        for i, etaval in enumerate(eta):
            krd = k0*MHC_kfunc(-eta_f[i], lmbda)
            kox = k0*MHC_kfunc(eta_f[i], lmbda)
            Rate[i] = ecd_extras[i]*(krd*c_lyte - kox*c_sld[i])
    else:
        krd = k0*MHC_kfunc(-eta_f, lmbda)
        kox = k0*MHC_kfunc(eta_f, lmbda)
        Rate = np.exp(-E_A/T + E_A/1) * ecd_extras*(krd*c_lyte - kox*c_sld)
    return Rate


def CIET(eta, c_sld, c_lyte, k0, E_A, T, act_R=None,
         act_lyte=None, lmbda=None, alpha=None):
    # See Fraggedakis et al. 2020
    eta_f = eta + T*np.log(c_lyte/c_sld)
    ecd_extras = (1-c_sld)/np.sqrt(4.0*np.pi*lmbda)
    if isinstance(eta, np.ndarray):
        Rate = np.empty(len(eta), dtype=object)
        for i, etaval in enumerate(eta):
            krd = k0*MHC_kfunc(-eta_f[i], lmbda)
            kox = k0*MHC_kfunc(eta_f[i], lmbda)
            Rate[i] = ecd_extras[i]*(krd*c_lyte - kox*c_sld[i])
    else:
        krd = k0*MHC_kfunc(-eta_f, lmbda)
        kox = k0*MHC_kfunc(eta_f, lmbda)
        Rate = np.exp(-E_A/T + E_A/1) * ecd_extras*(krd*c_lyte - kox*c_sld)
    return Rate
