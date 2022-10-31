import numpy as np


def dLiFePO4dc(self, y, ybar, muR_ref):
    """ Bai, Cogswell, Bazant 2011 """
    muRtheta = -self.eokT*3.422
    #muRhomog = self.reg_sln(y, self.get_trode_param("Omega_a"))
    #muRnonHomog = self.general_non_homog(y, ybar)
    #muR = muRhomog + muRnonHomog
    omega = self.get_trode_param("Omega_a")
    muR = (1-2*y*omega+2*y**2*omega)/(y-y**2)
    actR = np.exp(muR/self.T)
   # muR += muRtheta + muR_ref
    muR = muR * 0
    return muR, actR
