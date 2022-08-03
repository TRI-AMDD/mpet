from mpet.config import constants


def SolidSol_fitted(self, y, ybar, muR_ref, ISfuncs=None):
    T = self.T
    x = y
    OCV_ref = (4.22882816e+00 +
               - 7.46680452e-01*x - 6.46487121e-01*x**2
               + 1.98317348e+02*x**3 - 3.38658221e+03*x**4
               + 2.55770620e+04*x**5 - 1.10532458e+05*x**6
               + 3.00933346e+05*x**7 - 5.36070100e+05*x**8
               + 6.24445330e+05*x**9 - 4.52755354e+05*x**10
               + 1.70441348e+05*x**11 - 2.11992531e+03*x**12
               - 2.26426671e+04*x**13 + 5.91243316e+03*x**14)
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
    OCV = OCV_ref + dUdT*(T-1)*constants.T_ref
    muR = self.get_muR_from_OCV(OCV, muR_ref)
    actR = None
    return muR, actR
