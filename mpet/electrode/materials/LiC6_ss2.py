from mpet.props_am import step_down


def LiC6_ss2(self, y, ybar, muR_ref, ISfuncs=None):
    """ Bernardi and Go 2011 """
    p1, p2, p3, p4 = (0.085, 0.120, 0.210, 3.5)
    sfac = 0.3
    OCV = (p1*step_down(y, 1., sfac*0.02)
           + (p2 - p1)*step_down(y, 0.5, 0.005)
           + (p3 - p2)*step_down(y, 0.1944, sfac*0.03571)
           + (p4 - p3)*step_down(y, 0., sfac*0.08333))
    muR = self.get_muR_from_OCV(OCV, muR_ref)
    actR = None
    return muR, actR
