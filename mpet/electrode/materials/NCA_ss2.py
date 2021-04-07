import numpy as np


def NCA_ss2(self, y, ybar, muR_ref, ISfuncs=None):
    """
    Li_q Ni(0.8)Co(0.15)Al(0.05)O2
    as a function of y. Here, y actually represents a practical
    utilization of 70% of the material, so the material is "empty"
    (y=0) when q=0.3 and full (y=1) when q=1.
    This function was obtained from a fit by Raymond B. Smith
    of Samsung data of a LiC6-NCA cell discharged at C/100.
    """
    OCV = (-self.kToe*np.log(y/(1-y))
           + 4.12178 - 0.2338*y - 1.24566*y**2 + 1.16769*y**3
           - 0.20745*y**4)
    muR = self.get_muR_from_OCV(OCV, muR_ref)
    actR = None
    return muR, actR
