import numpy as np


def ACR_twoConc(self, y, ybar, muR_ref, ISfuncs=None):
    """ New test material """
    muRtheta1 = -self.eokT*4.05
    muRtheta2 = -self.eokT*3.5
    if ISfuncs is None:
        ISfuncs1, ISfuncs2 = None, None
    else:
        ISfuncs1, ISfuncs2 = ISfuncs
    y1, y2 = y
    Omga = self.get_trode_param('Omega_a')
    Omgb = self.get_trode_param('Omega_b')
    Omgc = self.get_trode_param('Omega_c')
    muR1 = self.reg_sln(y1, Omga, ISfuncs1)
    muR2 = self.reg_sln(y2, Omga, ISfuncs2)
    muR1nonHomog, muR2nonHomog = self.general_non_homog(y, ybar)
    # muR1 += (Omgb*(1-2*y1)*y2*(1-y2) + Omgc*(1-2*y2) - Omgc)
    # muR2 += (Omgb*(1-2*y2)*y1*(1-y1) + Omgc*(1-2*y1) - Omgc)
    muR1 += (Omgb*(6*y1**2 - 6*y1 + 1)
             + Omgc*((1-2*y1)**3-4*y1*(1-y1)*(1-2*y1)))
    muR2 += (Omgb*(6*y2**2 - 6*y2 + 1)
             + Omgc*((1-2*y2)**3-4*y2*(1-y2)*(1-2*y2)))
    muR1 += muR1nonHomog
    muR2 += muR2nonHomog
    actR1 = np.exp(muR1/self.T)
    actR2 = np.exp(muR2/self.T)
    muR1 += muRtheta1 + muR_ref
    muR2 += muRtheta2 + muR_ref
    return (muR1, muR2), (actR1, actR2)
