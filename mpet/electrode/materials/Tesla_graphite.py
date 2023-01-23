import numpy as np


def Tesla_graphite(self, y, ybar, muR_ref):
    """ Berliner et al., 2022.
    OCV for graphite [V]
    theta_p is the solid particle surface concentration normalized by its maximum concentration [-]
    """
    a0 = -48.9921992984694
    a = np.array([29.9816001180044, 161.854109570929, -0.283281555638378,
                 - 47.7685802868867, -65.0631963216785]).reshape([1,-1])
    b = np.array([0.005700461982903098, -0.1056830819588037, 0.044467320399373095,
                 - 18.947769999614668, 0.0022683366694012178]).reshape([1,-1])
    c = np.array([-0.050928145838337484, 0.09687316296868148, 0.04235223640014242,
                 7.040771011524739, 0.0011604439514018858]).reshape([1,-1])
    y = y.reshape([-1,1])

    OCV = a0 + np.squeeze(a[0,0]*np.exp((y - b[0,0])/c[0,0])) + \
        np.sum(a[0,1:]*np.tanh((y - b[0,1:])/c[0,1:]), axis=1)
    muR = self.get_muR_from_OCV(OCV, muR_ref)
    actR = None
    return muR, actR
