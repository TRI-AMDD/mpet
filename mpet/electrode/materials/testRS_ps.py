import numpy as np


def testRS_ps(self, y, ybar, muR_ref, ISfuncs=None):
    muRtheta = -self.eokT*2.
    muRhomog = self.reg_sln(y, self.ndD["Omga"], ISfuncs)
    muRnonHomog = self.general_non_homog(y, ybar)
    muR = muRhomog + muRnonHomog
    actR = np.exp(muR/self.T)
    muR += muRtheta + muR_ref
    return muR, actR
