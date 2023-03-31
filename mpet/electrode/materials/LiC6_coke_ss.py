import numpy as np


def LiC6_coke_ss(self, y, ybar, T, muR_ref):
    """ Doyle, Newman, 1996 """
    OCV = (-0.16 + 1.32*np.exp(-3.0*y) + 10.*np.exp(-2000.*y))
    muR = self.get_muR_from_OCV(OCV, muR_ref)
    actR = None
    return muR, actR
