import numpy as np


def LiFeMnPO4(self, y, ybar, muR_ref, ISfuncs=None):
    """ LFMP chemical potential """
    muRtheta1 = -self.eokT*4.09
    muRtheta2 = -self.eokT*3.422
    stoich_1 = self.get_trode_param("stoich_1")
    stoich_2 = 1-stoich_1
    if ISfuncs is None:
        ISfuncs1, ISfuncs2 = None, None
    else:
        ISfuncs1, ISfuncs2 = ISfuncs
    y1, y2 = y
    Omga = self.get_trode_param('Omega_a')
    Omgb = self.get_trode_param('Omega_b')
    Omgc = (Omga+Omgb)*0.3
    # 2*ln(c/1-c) is used to account for electrons config entropy
    muR1homog = 2*self.reg_sln(y1, (Omga/2)*stoich_1, ISfuncs1)
    muR2homog = 2*self.reg_sln(y2, (Omgb/2)*stoich_2, ISfuncs2)
    muR1nonHomog, muR2nonHomog = self.general_non_homog(y, ybar)
    muR1 = muR1homog + muR1nonHomog
    muR2 = muR2homog + muR2nonHomog
    # interaction between the two phases
    muR1 += stoich_2*Omgc*(1-2*y2)
    muR2 += stoich_1*Omgc*(1-2*y1)
    actR1 = np.exp(muR1/self.T)
    actR2 = np.exp(muR2/self.T) 
    muR1 += muRtheta1 + muR_ref
    muR2 += muRtheta2 + muR_ref
    return (muR1, muR2), (actR1, actR2)
