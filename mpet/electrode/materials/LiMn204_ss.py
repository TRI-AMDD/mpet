import numpy as np


def LiMn2O4_ss(self, y, ybar, muR_ref, ISfuncs=None):
    """ Doyle, Newman, 1996 """
    # OCV in V vs Li/Li+
    OCV = (4.19829 + 0.0565661*np.tanh(-14.5546*y + 8.60942)
           - 0.0275479*(1/((0.998432 - y)**(0.492465)) - 1.90111)
           - 0.157123*np.exp(-0.04738*y**8)
           + 0.810239*np.exp(-40*(y - 0.133875)))
    muR = self.get_muR_from_OCV(OCV, muR_ref)
    actR = None
    return muR, actR
