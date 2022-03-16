import numpy as np


def testRS(self, y, ybar, T, muR_ref, ISfuncs=None):
    muRtheta = 0.
    muR = self.reg_sln(y, T, self.get_trode_param("Omega_a"), ISfuncs)
    actR = np.exp(muR/T)
    muR += muRtheta + muR_ref
    return muR, actR
