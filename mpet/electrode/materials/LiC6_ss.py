import numpy as np


def LiC6_ss(self, y, ybar, muR_ref, ISfuncs=None):
    """ Safari, Delacourt 2011 """
    OCV = (0.6379 + 0.5416*np.exp(-305.5309*y)
           + 0.044*np.tanh(-(y - 0.1958)/0.1088)
           - 0.1978*np.tanh((y - 1.0571)/0.0854)
           - 0.6875*np.tanh((y + 0.0117)/0.0529)
           - 0.0175*np.tanh((y - 0.5692)/0.0875))
    muR = self.get_muR_from_OCV(OCV, muR_ref)
    actR = None
    return muR, actR
