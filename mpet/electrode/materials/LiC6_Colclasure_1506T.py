import numpy as np


def LiC6_Colclasure_1506T(self, y, ybar, muR_ref, ISfuncs=None):
    """Taken from Colclasure 2020, but only for lithiating. The range of\
    confidence is 0.01 to 0.97 for voltages of ~0.6V and \
    0.0435V. Don’t calculate U2 for intercalation fractions below 0.8, \
    it will cause numeric issues just use U1."""
    x = y
    U1 = \
        - 1.059423355572770E-02*np.tanh((x - 1.453708425609560E-02) / 9.089868397988610E-05) \
        + 2.443615203087110E-02*np.tanh((x - 5.464261369950400E-01) / 6.270508166379020E-01) \
        - 1.637520788053810E-02*np.tanh((x - 5.639025014475490E-01) / 7.053886409518520E-02) \
        - 6.542365622896410E-02*np.tanh((x - 5.960370524233590E-01) / 1.409966536648620E+00) \
        - 4.173226059293490E-02*np.tanh((x - 1.787670587868640E-01) / 7.693844911793470E-02) \
        - 4.792178163846890E-01*np.tanh((x + 3.845707852011820E-03) / 4.112633446959460E-02) \
        + 6.594735004847470E-01 \
        - 4.364293924074990E-02*np.tanh((x - 9.449231893318330E-02) / -2.046776012570780E-02) \
        - 8.241166396760410E-02*np.tanh((x - 7.746685789572230E-02) / 3.593817905677970E-02)
    U2 = \
        -1.731504647676420E+02*np.power(x, 8.0) + 8.252008712749000E+01*np.power(x, 7.0) \
        + 1.233160814852810E+02*np.power(x, 6.0) + 5.913206621637760E+01*np.power(x, 5.0) \
        + 3.322960033709470E+01*np.power(x, 4.0) + 3.437968012320620E+00*np.power(x, 3.0) \
        - 6.906367679257650E+01*np.power(x, 2.0) - \
        1.228217254296760E+01*x - 5.037944982759270E+01
    U = U1 + (U2 - U1) / (1.0 + np.exp(-1.0e2*(x - 1.02956203215198)))
    muR = self.get_muR_from_OCV(U, muR_ref)
    actR = None
    return muR, actR
