import numpy as np


def LiC6_1param(self, y, ybar, muR_ref, ISfuncs=None):
    muRtheta = -self.eokT*0.12
    config = self.config
    muRhomog = self.graphite_1param_homog_3(
        y, config["Omga"], config["Omgb"], ISfuncs)
    muRnonHomog = self.general_non_homog(y, ybar)
    muR = muRhomog + muRnonHomog
    actR = np.exp(muR/self.T)
    muR += muRtheta + muR_ref
    return muR, actR
