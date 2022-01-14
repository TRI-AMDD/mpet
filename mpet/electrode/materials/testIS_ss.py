import numpy as np


def testIS_ss(self, y, ybar, muR_ref, ISfuncs=None):
    """Ideal solution material for testing."""
    OCV = -self.kToe*np.log(y/(1-y))
    muR = self.get_muR_from_OCV(OCV, muR_ref)
    actR = None
    return muR, actR
