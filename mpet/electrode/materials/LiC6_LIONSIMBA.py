import numpy as np


def LiC6_LIONSIMBA(self, y, ybar, muR_ref, ISfuncs=None):
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
    OCV_ref = r1 + r2*y + r3*y**0.5 - r4 * \
        y**(-1) + r5*y**(-1.5) + r6*np.exp(0.9-15*y) - r7*np.exp(0.4465*y-0.4108)
    k1 = 0.001
    k2 = 0.005269056
    k3 = 3.299265709
    k4 = -91.79325798
    k5 = 1004.911008
    k6 = -5812.278127
    k7 = 19329.7549
    k8 = -37147.8947
    k9 = 38379.18127
    k10 = -16515.05308
    k11 = 1
    k12 = -48.09287227
    k13 = 1017.234804
    k14 = -10481.80419
    k15 = 59431.3
    k16 = -195881.6488
    k17 = 374577.3152
    k18 = -385821.1607
    k19 = 165705.8597
    dUdT = k1*(k2+k3*y+k4*y**2+k5*y**3+k6*y**4+k7*y**5+k8*y**6+k9*y**7+k10*y**8) / \
        (k11+k12*y+k13*y**2+k14*y**3+k15*y**4+k16*y**5+k17*y**6+k18*y**7+k19*y**8)
    OCV = OCV_ref + dUdT*(T-1)*Tref
    muR = self.get_muR_from_OCV(OCV, muR_ref)
    actR = None
    return muR, actR
