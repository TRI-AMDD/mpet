import numpy as np


def dLiC6_LIONSIMBAdc(self, y, ybar, muR_ref):
    """ Torchio et al, 2016. """
    T = self.T
    Tref = 298
    r1 = 0.7222
    r2 = 0.1387
    r3 = 0.029
    r4 = 0.0172
    r5 = 0.0019
    r6 = 0.2808
    r7 = 0.7984
    OCV_ref = r2 - 15*r6*np.exp(9/10 - 15*y) - (893*r7*np.exp((893*y)/2000 - 1027/2500)
                                                )/2000 + r3/(2*y**(1/2)) + r4/y**2 - (3*r5)/(2*y**(5/2))
    muR = self.get_muR_from_OCV(OCV_ref, 0)
    actR = None
    return muR, actR
