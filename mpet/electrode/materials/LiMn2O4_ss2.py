import numpy as np


def LiMn2O4_ss2(self, y, ybar, muR_ref, ISfuncs=None):
    """ Fuller, Doyle, Newman, 1994 """
    # OCV in V vs Li/Li+
    OCV = (4.06279 + 0.0677504*np.tanh(-21.8502*y + 12.8268)
           - 0.105734*(1/((1.00167 - y)**(0.379571)) - 1.575994)
           - 0.045*np.exp(-71.69*y**8)
           + 0.01*np.exp(-200*(y - 0.19)))
    muR = self.get_muR_from_OCV(OCV, muR_ref)
    actR = None
    return muR, actR
