import numpy as np


def LiC6_coke_ss2(self, y, ybar, muR_ref, ISfuncs=None):
    """ Fuller, Doyle, Newman, 1994 """
    c1 = -0.132056
    c2 = 1.40854
    c3 = -3.52312
    OCV = c1 + c2*np.exp(c3*y)
    muR = self.get_muR_from_OCV(OCV, muR_ref)
    actR = None
    return muR, actR
