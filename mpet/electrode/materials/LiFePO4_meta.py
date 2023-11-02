import numpy as np


def LiFePO4_meta(self, y, ybar, T, muR_ref):
    muRtheta = -self.eokT*3.422
    k = 4
    h = 40
    b = 0.1
    a = -0.03
    phi = -40
    muRhomog = self.reg_sln(y, T, self.get_trode_param("Omega_a"))
    muRhomog += k*phi*(1-2*y)*(y*(1-y))**(k-1)
    muRhomog += h*a*np.pi*np.cos(np.pi*(y-b))*np.sin(np.pi*(y-b))**(h-1)
    muRnonHomog = self.general_non_homog(y, ybar)
    muR = muRhomog + muRnonHomog
    actR = np.exp(muR/T)
    muR += muRtheta + muR_ref
    return muR, actR
