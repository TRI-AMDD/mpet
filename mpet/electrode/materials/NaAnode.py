import numpy as np

#     muRtheta = -self.eokT*1.3356
#     OCV = (-13.522*y**6 + 42.451*y**5 - 46.305*y**4
#            + 15.214*y**3 + 7.0034*y**2 - 6.0744*y)


def NaAnode(self, y, ybar, T, muR_ref):

    muRtheta = -self.eokT*1.435080161287005
    muRtheta = -self.eokT*1.45

    OCV = (-7.6991314610297525*y
           + 16.406426373526156*y**2 + -13.120937152671123*y**3
           + -8.140487082887885*y**4 + 22.754536573994052*y**5
           + -11.76391247146248*y**6)
    
#     muRtheta = -self.eokT*1.2356
#     OCV = (-13.522*y**6 + 42.451*y**5 - 46.305*y**4
#            + 15.214*y**3 + 7.0034*y**2 - 6.0744*y)

    muRhomog = -self.eokT*OCV
    muRnonHomog = self.general_non_homog(y, ybar)
    muR = muRhomog + muRnonHomog
    actR = np.exp(muR/T)
    muR += muRtheta + muR_ref

    return muR, actR
