import numpy as np


def LiFePO4_rs_red(self, y, ybar, T, muR_ref):
    muRtheta = -self.eokT*3.422
    muRhomog = self.reg_sol_red(y, T, self.get_trode_param("Omega_a"),muRtheta)
    muRnonHomog = self.general_non_homog(y, ybar)
    muR = muRhomog + muRnonHomog
    actR = np.exp(muR/T)
    muR += muR_ref
    return muR, actR
