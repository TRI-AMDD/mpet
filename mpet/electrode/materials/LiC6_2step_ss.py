import numpy as np
from mpet.props_am import step_up, step_down


def LiC6_2step_ss(self, y, ybar, muR_ref, ISfuncs=None):
    """
    Fit function to the OCV predicted by the phase separating
    2-variable graphite model (LiC6 function in this class).
    """
    Vstd = 0.1196
    Vstep = 0.0351733976
    edgeLen = 0.024
    lEdge = edgeLen
    rEdge = 1 - edgeLen
    width = 1e-4
    vshift = 1e-2
    lSide = -((np.log(y/(1-y)) - np.log(lEdge/(1-lEdge)) - vshift)*step_down(y, lEdge, width))
    rSide = -((np.log(y/(1-y)) - np.log(rEdge/(1-rEdge)) + vshift)*step_up(y, rEdge, width))
    OCV = (
        Vstd
        + Vstep*(step_down(y, 0.5, 0.013) - 1)
        + self.kToe*(lSide + rSide)
        )
    muR = self.get_muR_from_OCV(OCV, muR_ref)
    actR = None
    return muR, actR
