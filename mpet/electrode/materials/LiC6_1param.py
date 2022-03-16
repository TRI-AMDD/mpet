import numpy as np


def LiC6_1param(self, y, ybar, T, muR_ref, ISfuncs=None):
    muRtheta = -self.eokT*0.12
    muRhomog = self.graphite_1param_homog_3(
        y, T, self.get_trode_param("Omega_a"), self.get_trode_param("Omega_b"), ISfuncs)
    muRnonHomog = self.general_non_homog(y, ybar)
    muR = muRhomog + muRnonHomog
    actR = np.exp(muR/T)
    muR += muRtheta + muR_ref
    return muR, actR
