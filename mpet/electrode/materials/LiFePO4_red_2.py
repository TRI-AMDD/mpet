import numpy as np
from mpet.props_am import step_up


def LiFePO4_red_2(self, y, ybar, T, muR_ref, ISfuncs=None):
    """ Huada: to capture tails """
    # muRtheta = -self.eokT*3.435
    # muRhomog = -15.0*np.exp(-y/0.005) + (1.5-2.0*y)*step_up(y,0.03,0.02) \
    #             + 20.0*step_up(y, 1.0, 0.02)
    muRtheta = -self.eokT*3.405
    # muRhomog = -15.0*np.exp(-(y-0.02)/0.02) + (1.5-4.3*y)*step_up(y,0.03,0.05) \
    #             + 20.0*step_up(y, 1.0, 0.07)
    muRhomog = -15.0*np.exp(-(y-0.02)/0.02) + (1.5-4.3*y)*step_up(y,0.03,0.05) \
                + 20.0*step_up(y, 1.0, 0.06)
    muRnonHomog = self.general_non_homog(y, ybar)
    muR = muRhomog + muRnonHomog
    actR = np.exp(muR/T)
    muR += muRtheta + muR_ref
    return muR, actR