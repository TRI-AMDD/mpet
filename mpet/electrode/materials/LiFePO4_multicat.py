import numpy as np


def LiFePO4_multicat(self, y, ybar, T, muR_ref):
    """ Ombrini 2023 """
    muRtheta1 = -self.eokT*0.1
    # muRtheta2 = -self.eokT*3.422
    muRtheta2 = -self.eokT*0.3
    y1, y2 = y
    Omga = self.get_trode_param('Omega_a')
    Omgb = self.get_trode_param('Omega_b')
    Omgc = self.get_trode_param('Omega_c')
    # ln(c/1-c) is used to account for electrons config entropy
    muR1homog = np.log(y1/(1-y2-y1)) + Omga*(1-y2-2*y1)
    muR2homog = np.log(y2/(1-y2-y1)) + Omgb*(1-y1-2*y2)
    muR1nonHomog, muR2nonHomog = self.general_non_homog(y, ybar)
    muR1 = muR1homog + muR1nonHomog
    muR2 = muR2homog + muR2nonHomog
    # interaction between the two phases
    muR1 += - Omgb*y2 + Omgc*y2
    muR2 += - Omga*y1 + Omgc*y1
    actR1 = np.exp(muR1/T)
    actR2 = np.exp(muR2/T)
    muR1 += muRtheta1 + muR_ref
    muR2 += muRtheta2 + muR_ref
    return (muR1, muR2), (actR1, actR2)
