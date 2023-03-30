import numpy as np


def LiFePO4_meta(self, y, ybar, T, muR_ref):
    """ Bai, Cogswell, Bazant 2011 """
    muRtheta = -self.eokT*3.422
    muRhomog = 2*self.reg_sln(y, T, self.get_trode_param("Omega_a")/2)
    b = 0.3
    pi = 3.14159
    a = -0.11
    phi = -700
    muR_meta_1 = a*20*pi*np.cos(pi*y-b)*(np.sin(pi*y-b))**19
    muR_meta_2 = phi*6*((y*(1-y))**5)*(1-2*y)
    # a = -0.1
    # muR_meta = a*4*pi*np.cos(pi*y-b)*(np.sin(pi*y-b))**3
    muRnonHomog = self.general_non_homog(y, ybar)
    muR = muRhomog + muRnonHomog + muR_meta_1 + muR_meta_2
    actR = np.exp(muR/T)
    muR += muRtheta + muR_ref
    return muR, actR
