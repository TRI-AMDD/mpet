import numpy as np


def SolidSold_fromEnt(self, y, ybar, muR_ref, ISfuncs=None):
    muRtheta = -self.eokT * 3.88
    muRhomog = self.reg_sln(y, self.get_trode_param("Omega_a"), ISfuncs)
    muRnonHomog = self.general_non_homog(y, ybar)
    muR = muRhomog + muRnonHomog
    # Omgb = self.get_trode_param("Omega_b")
    # muR += Omgb*((1-2*y)**2 - 2*y*(1-y))
    actR = np.exp(muR / self.T)
    muR += muRtheta + muR_ref
    return muR, actR
