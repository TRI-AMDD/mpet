from mpet.props_am import step_down, step_up


def testRS_ss(self, y, ybar, muR_ref, ISfuncs=None):
    """
    Regular solution material which phase separates at binodal points,
    for modeling as a solid solution. For testing.
    """
    # Based Omg = 3*k*T_ref
    yL = 0.07072018
    yR = 0.92927982
    OCV_rs = -self.kToe*self.reg_sln(y, self.get_trode_param("Omega_a"), ISfuncs)
    width = 0.005
    OCV = OCV_rs*step_down(y, yL, width) + OCV_rs*step_up(y, yR, width) + 2
    muR = self.get_muR_from_OCV(OCV, muR_ref)
    actR = None
    return muR, actR
