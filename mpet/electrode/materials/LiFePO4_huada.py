import numpy as np
from mpet.props_am import step_down, step_up


def LiFePO4_red(self, y, ybar, T, muR_ref):
    """ Huada: to capture tails """
    muRtheta = -self.eokT*3.42
    muRhomog = -15.0*np.exp(-y/0.005) + (1.5-2.0*y)*step_up(y,0.03,0.02) \
                + 20.0*step_up(y, 1.0, 0.02)
    muRnonHomog = self.general_non_homog(y, ybar)
    muR = muRhomog + muRnonHomog
    actR = np.exp(muR/T)
    muR += muRtheta + muR_ref
    return muR, actR