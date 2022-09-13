import numpy as np


def testRS(self, y, ybar, muR_ref):
    muRtheta = 0.
    muR = self.reg_sln(y, self.get_trode_param("Omega_a"))
    actR = np.exp(muR/self.T)
    muR += muRtheta + muR_ref
    return muR, actR
