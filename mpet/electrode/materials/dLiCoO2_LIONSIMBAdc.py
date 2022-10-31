def dLiCoO2_LIONSIMBAdc(self, y, ybar, muR_ref):
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
    OCV_ref = (2*y*(5*r12*y**8 - 4*r11*y**6 + 3*r10*y**4 - 2*r9*y**2 + r8)*(- r6*y**10 + r5*y**8 - r4*y**6 + r3*y**4 - r2*y**2 + r1))/(- r12*y**10 + r11*y**8 - \
               r10*y**6 + r9*y**4 - r8*y**2 + r7)**2 - (2*y*(5*r6*y**8 - 4*r5*y**6 + 3*r4*y**4 - 2*r3*y**2 + r2))/(- r12*y**10 + r11*y**8 - r10*y**6 + r9*y**4 - r8*y**2 + r7)
    muR = self.get_muR_from_OCV(OCV_ref, 0)
    actR = None
    return muR, actR
