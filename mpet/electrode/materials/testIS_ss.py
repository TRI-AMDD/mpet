import numpy as np


def testIS_ss(self, y, ybar, T, muR_ref):
    """Ideal solution material for testing."""
    OCV = -self.kToe*np.log(y/(1-y))*T
    muR = self.get_muR_from_OCV(OCV, muR_ref)
    actR = None
    return muR, actR
