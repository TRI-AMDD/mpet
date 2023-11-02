import numpy as np


def LiFeMnPO4_red(self, y, ybar, T, muR_ref):
    """ Ombrini 2023 """
    muRtheta1 = -self.eokT*4.09
    muRtheta2 = -self.eokT*3.422
    stoich_1 = self.get_trode_param("stoich_1")
    Omga = self.get_trode_param('Omega_a')
    Omgb = self.get_trode_param('Omega_b')
    muRhomog = self.LMFP_red(y, T, Omga, Omgb, muRtheta1, muRtheta2, stoich_1)
    muRnonHomog = self.general_non_homog(y, ybar)
    muR = muRhomog + muRnonHomog
    actR = np.exp(muR/T)
    muR += muR_ref
    return muR, actR
