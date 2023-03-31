import numpy as np


def testRS_ps(self, y, ybar, T, muR_ref, ISfuncs=None):
    muRtheta = -self.eokT*2.
    muRhomog = self.reg_sln(y, T, self.get_trode_param("Omega_a"), ISfuncs)
    muRnonHomog = self.general_non_homog(y, ybar)
    muR = muRhomog + muRnonHomog
    actR = np.exp(muR/T)
    muR += muRtheta + muR_ref
    return muR, actR
