import numpy as np


def LiFePO4_meta(self, y, ybar, T, muR_ref):
    """ Bai, Cogswell, Bazant 2011 """
    muRtheta = -self.eokT*3.422
    muRhomog = 2*self.reg_sln(y, T, self.get_trode_param("Omega_a")/2)
    a = -0.2
    b = 0.15
    pi = 3.14159
    muR_meta = a*6*pi*np.cos(pi*y-b)*(np.sin(pi*y-b))**5
    muRnonHomog = self.general_non_homog(y, ybar)
    muR = muRhomog + muRnonHomog + muR_meta
    actR = np.exp(muR/T)
    muR += muRtheta + muR_ref
    return muR, actR
