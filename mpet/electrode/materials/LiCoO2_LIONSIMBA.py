def LiCoO2_LIONSIMBA(self, y, ybar, muR_ref, ISfuncs=None):
    """ Torchio et al, 2016. """
    T = self.T
    Tref = 298
    r1 = 4.656
    r2 = 88.669
    r3 = 401.119
    r4 = 342.909
    r5 = 462.471
    r6 = 433.434
    r7 = 1
    r8 = 18.933
    r9 = 79.532
    r10 = 37.311
    r11 = 73.083
    r12 = 95.96
    OCV_ref = (-r1 + r2*y**2 - r3*y**4 + r4*y**6 - r5*y**8 + r6 * y**10) / \
        (-r7 + r8*y**2 - r9*y**4 + r10*y**6 - r11*y**8 + r12*y**10)
    k1 = -0.001
    k2 = 0.199521039
    k3 = -0.928373822
    k4 = 1.364550689000003
    k5 = -0.6115448939999998
    k6 = 1
    k7 = -5.661479886999997
    k8 = 11.47636191
    k9 = -9.82431213599998
    k10 = 3.048755063
    dUdT = k1*(k2+k3*y+k4*y**2+k5*y**3)/(k6+k7*y+k8*y**2+k9*y**3+k10*y**4)
    OCV = OCV_ref + dUdT*(T-1)*Tref
    muR = self.get_muR_from_OCV(OCV, muR_ref)
    actR = None
    return muR, actR
