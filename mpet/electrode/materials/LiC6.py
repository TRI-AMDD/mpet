import numpy as np


def LiC6(self, y, ybar, T, muR_ref, ISfuncs=(None, None)):
    """ Ferguson and Bazant 2014 """
    muRtheta = -self.eokT*0.12
    muR1homog, muR2homog = self.graphite_2param_homog(
        y, T, self.get_trode_param("Omega_a"), self.get_trode_param("Omega_b"),
        self.get_trode_param("Omega_c"), self.get_trode_param("EvdW"), ISfuncs)
    muR1nonHomog, muR2nonHomog = self.general_non_homog(y, ybar)
    muR1 = muR1homog + muR1nonHomog
    muR2 = muR2homog + muR2nonHomog
    actR1 = np.exp(muR1/T)
    actR2 = np.exp(muR2/T)
    muR1 += muRtheta + muR_ref
    muR2 += muRtheta + muR_ref
    return (muR1, muR2), (actR1, actR2)
