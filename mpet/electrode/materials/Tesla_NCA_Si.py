import numpy as np


def Tesla_NCA_Si(self, y, ybar, muR_ref):
    """ Berliner et al., 2022.
    OCV for graphite [V]
    theta_p is the solid particle surface concentration normalized by its maximum concentration [-]
    """
    a = np.array([0.145584993881910,2.526321858618340,172.0810484337340,1.007518156438100,
                 1.349501707184530,0.420519124096827,2.635800979146210,
                 3.284611867463240]).reshape([1,-1])
    b = np.array([0.7961299985689542, 0.2953029849791878, -1.3438627370872127, 0.6463272973815986,
                 0.7378056244779166, 0.948857021183584, 0.5372357238527894,
                 0.8922020984716097]).reshape([1,-1])
    c = np.array([0.060350976183950786, 0.20193410562543265, 0.7371221766768185,
                 0.10337785458522612, 0.09513470475980132, 0.0422930728072207, 0.1757549310633964,
                 0.1413934223088055]).reshape([1,-1])
    y = y.reshape([-1,1])

    OCV = np.sum(a*np.exp(-((y - b)/c)**2), axis=1)
    muR = self.get_muR_from_OCV(OCV, muR_ref)
    actR = None
    return muR, actR
