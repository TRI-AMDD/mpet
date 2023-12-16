import numpy as np

def mu(x, y, Lxs, Lys, Lxy, Lxyv):
    mu = np.log(x / (1-y-x))
    kx = 0
    for Lx in Lxs:
        mu += Lx * ((1-y-2*x)**kx) * ((1-y-2*x) - 2 * kx * (x * (1-y-x)/(1-y-2*x)))
        kx += 1
    ky = 0
    for Ly in Lys:
        mu += - Ly * y * ((1-x-2*y)**ky) * (1 + ky*(1-x-y)/(1-x-2*y))
        ky += 1
    mu += Lxy*y
    mu += Lxyv*y*(1-y-2*x)
    return mu


def LiNaFePO4(self, y, ybar, T, muR_ref):
    """ Ombrini 2024 """
    muRtheta1 = -self.eokT*0.3
    muRtheta2 = -self.eokT*0.1
    y1, y2 = y
    Omga = self.get_trode_param('Omega_a')
    Lnali = self.get_trode_param('Omega_b')
    Lnalivac = self.get_trode_param('Omega_c')
    L_Li = [Omga]
    L_Na = [0.94502646, 8.02136736, 5.35420982, -15.21264346, -4.0081790, 7.62295359]

    muLihom = mu(y1, y2, L_Li, L_Na, Lnali, Lnalivac)
    muNahom = mu(y2, y1, L_Na, L_Li, Lnali, Lnalivac)

    muR1nonHomog, muR2nonHomog = self.general_non_homog(y, ybar)
    muR1 = muLihom + muR1nonHomog
    muR2 = muNahom + muR2nonHomog
    # interaction between the two phases
    actR1 = np.exp(muR1/T)
    actR2 = np.exp(muR2/T)
    muR1 += muRtheta1 + muR_ref
    muR2 += muRtheta2 + muR_ref
    return (muR1, muR2), (actR1, actR2)


