import numpy as np


def NCA_ss1(self, y, ybar, muR_ref, ISfuncs=None):
    """
    This function was obtained from Dan Cogswell's fit of Samsung
    data.
    """
    OCV = (3.86 + 1.67*y - 9.52*y**2 + 15.04*y**3 - 7.95*y**4
           - 0.06*np.log(y/(1-y)))
    muR = self.get_muR_from_OCV(OCV, muR_ref)
    actR = None
    return muR, actR
