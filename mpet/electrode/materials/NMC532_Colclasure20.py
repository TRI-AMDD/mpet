import numpy as np


def NMC532_Colclasure20(self, y, ybar, muR_ref, ISfuncs=None):
    x = y
    OCV = (5.314735633000300E+00 +
           -3.640117692001490E+03*x**14.0 + 1.317657544484270E+04*x**13.0
           - 1.455742062291360E+04*x**12.0 - 1.571094264365090E+03*x**11.0
           + 1.265630978512400E+04*x**10.0 - 2.057808873526350E+03*x**9.0
           - 1.074374333186190E+04*x**8.0 + 8.698112755348720E+03*x**7.0
           - 8.297904604107030E+02*x**6.0 - 2.073765547574810E+03*x**5.0
           + 1.190223421193310E+03*x**4.0 - 2.724851668445780E+02*x**3.0
           + 2.723409218042130E+01*x**2.0 - 4.158276603609060E+00*x +
           -5.573191762723310E-04*np.exp(6.560240842659690E+00*x**4.148209275061330E+01)
           )
    muR = self.get_muR_from_OCV(OCV, muR_ref)
    actR = None
    return muR, actR
