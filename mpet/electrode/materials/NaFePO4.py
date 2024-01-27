import numpy as np

def mu(x, Lxs):
    mu = np.log(x / (1-x))
    kx = 0
    for Lx in Lxs:
        mu += Lx * ((1-2*x)**kx) * ((1-2*x) - 2 * kx * (x * (1-x)/(1-2*x)))
        kx += 1
    return mu


def NaFePO4(self, y, ybar, T, muR_ref):
    muRtheta = -self.eokT*0
    L_Na = [0.94502646, 8.02136736, 5.35420982, -15.21264346, -4.0081790, 7.62295359]
    muRhomog = mu(y, L_Na)
    muRnonHomog = self.general_non_homog(y, ybar)
    muR = muRhomog + muRnonHomog
    actR = np.exp(muR/T)
    muR += muRtheta + muR_ref
    return muR, actR
