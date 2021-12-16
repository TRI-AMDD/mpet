import numpy as np


def LiC6(self, y, ybar, muR_ref, ISfuncs=(None, None)):
    """ Ferguson and Bazant 2014 """
    muRtheta = -self.eokT*0.12
    config = self.config
    muR1homog, muR2homog = self.graphite_2param_homog(
        y, config["Omga"], config["Omgb"], config["Omgc"], config["EvdW"], ISfuncs)
    muR1nonHomog, muR2nonHomog = self.general_non_homog(y, ybar)
    muR1 = muR1homog + muR1nonHomog
    muR2 = muR2homog + muR2nonHomog
    actR1 = np.exp(muR1/self.T)
    actR2 = np.exp(muR2/self.T)
    muR1 += muRtheta + muR_ref
    muR2 += muRtheta + muR_ref
    return (muR1, muR2), (actR1, actR2)
