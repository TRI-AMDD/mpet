import numpy as np


def LiFePO4_aniso(self, y, ybar, T, muR_ref):
    """ Ombrini 2023 """
    muRtheta1 = -self.eokT*3.422
    muRtheta2 = -self.eokT*3.422
    defects = np.random.rand(np.size(y[0]))*0.1
    # stoich_1 = self.get_trode_param("stoich_1")
    # stoich_2 = 1-stoich_1
    y1, y2 = y
    Omga = self.get_trode_param('Omega_a')
    Omgb = self.get_trode_param('Omega_b')
    Omgc = Omga - Omgb
    # 2*ln(c/1-c) is used to account for electrons config entropy
    # muR1homog = self.reg_sln(y1, T, Omgc)
    muR1homog = T*np.log(y1/(1-defects-y1)) + Omgc*(1-defects-2*y1)
    # muR2homog = self.reg_sln(y2, T, Omgc)
    muR2homog = T*np.log(y2/(1-defects-y2)) + Omgc*(1-defects-2*y2)
    muR1nonHomog, muR2nonHomog = self.general_non_homog(y, ybar)
    muR1 = muR1homog + muR1nonHomog
    muR2 = muR2homog + muR2nonHomog
    # interaction between the two phases
    muR1 += Omgb*(1-defects-2*y2)
    muR2 += Omgb*(1-defects-2*y1)
    actR1 = np.exp(muR1/T)
    actR2 = np.exp(muR2/T)
    muR1 += muRtheta1 + muR_ref
    muR2 += muRtheta2 + muR_ref
    return (muR1, muR2), (actR1, actR2)
