import numpy as np


def LTO(self, y, ybar, muR_ref, ISfuncs=None):
    """
    Vasileiadis 2017
    """
    muRtheta = -self.eokT * 1.55
    muRhomog = self.reg_sln(y, self.get_trode_param("Omega_a"))
    muRnonHomog = self.general_non_homog(y, ybar)
    muR = muRhomog + muRnonHomog
    actR = np.exp(muR / self.T)
    muR += muRtheta + muR_ref
    return muR, actR
