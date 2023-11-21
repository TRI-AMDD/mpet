import numpy as np
from scipy.ndimage import gaussian_filter

def LiFePO4_def(self, y, ybar, T, muR_ref):
    # d = np.random.rand(np.size(y))*0.3
    mean = 0.02
    stddevs = 0.02
    var = stddevs**2
    mu = np.log(mean**2/np.sqrt(var+mean**2))
    sigma = np.sqrt(np.log(var/mean**2 + 1))
    d = np.random.lognormal(mu, sigma,np.size(y))
    # gaussian smoothing
    d = gaussian_filter(d, sigma=0.5)

    muRtheta = -self.eokT*3.422 
    muRhomog = T*np.log(y/(1-d-y))+self.get_trode_param("Omega_a")*(1-d-2*y)
    # muRhomog = self.reg_sln(y, T, self.get_trode_param("Omega_a"))
    muRnonHomog = self.general_non_homog(y, ybar)
    muR = muRhomog + muRnonHomog
    actR = np.exp(muR/T)
    muR += muRtheta + muR_ref
    return muR, actR
