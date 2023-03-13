import numpy as np


def NaAnode(self, y, ybar, T, muR_ref):
    v0 = 1.4356
    OCV = (-13.522*y**6 + 42.451*y**5 - 46.305*y**4
           + 15.214*y**3 + 7.0034*y**2 - 6.0744*y + v0)
    muR = self.get_muR_from_OCV(OCV, muR_ref)
    actR = np.exp(muR/T)
    return muR, actR
