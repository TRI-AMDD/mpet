import numpy as np

import mpetMaterials

class muRfuncs():
    """ This class defines functions which describe the chemical
    potential of materials.
    Each function describes a particular material. The
    chemical potential is directly related to the open circuit voltage
    (OCV) for solid solution materials.
    In each function, muR_ref is the offset (non-dimensional) chemical
    potential by which the returned value will be shifted (for
    numerical convenience).
    Each material function returns:
        muR -- chemical potential
        actR -- activity (if applicable, else None)
    """
    def __init__(self, T, ndD=None, **kwargs):
        """ ndD can be the full dictionary of nondimensional
        parameters for the electrode particles, as made for the
        simulations. Otherwise, parameters can be passed directly in
        as keyword arguments and must contain all of the needed
        parameters for the material of interest.
        E.g.
        For a regular solution material:
            muRfuncs(T, ndD)
        or
            muRfuncs(T, muRfunc="LiFePO4", Omga=3.4)
        For solid solution function based on fit OCV:
            muRfuncs(T, ndD)
        or
            muRfuncs(T, muRfunc="LiMn2O4_ss")
        """
        if ndD is None:
            ndD = kwargs
        self.ndD = ndD
        self.T = T # nondimensional
        k = 1.381e-23
        Tabs = 298
        e = 1.602e-19
        self.eokT = e/(k*Tabs)
        self.kToe = (k*Tabs)/e
        materialData = {}
        materialData['LiMn2O4_ss'] = self.LiMn2O4_ss
        materialData['LiMn2O4_ss2'] = self.LiMn2O4_ss2
        materialData['LiC6_coke_ss'] = self.LiC6_coke_ss
        materialData['LiC6_coke_ss2'] = self.LiC6_coke_ss2
        materialData['LiC6_ss'] = self.LiC6_ss
        materialData['LiC6_ss2'] = self.LiC6_ss2
        materialData['LiC6_2step_ss'] = self.LiC6_2step_ss
        materialData['Li_ss'] = self.Li_ss
        materialData['NCA_ss1'] = self.NCA_ss1
        materialData['NCA_ss2'] = self.NCA_ss2
        materialData['LiFePO4'] = self.LiFePO4
        materialData['LiC6'] = self.LiC6
        self.muRfunc = materialData[ndD["muRfunc"]]

    def get_muR_from_OCV(self, OCV, muR_ref):
        return -self.eokT*OCV + muR_ref
    def get_actR_None(self, y):
        try:
            actR = np.array(y.shape[0]*[None])
        except AttributeError:
            actR = None
        return actR

    ######
    # Solid solution functions
    # These are all obtained directly from fitting an OCV of the
    # material with a standard counter electrode.
    # They can all only return values at 298 K
    ######

    def LiMn2O4_ss(self, y, ybar, muR_ref, ISfuncs=None):
        """ Doyle, Newman, 1996 """
        # OCV in V vs Li/Li+
        OCV = (4.19829 + 0.0565661*np.tanh(-14.5546*y + 8.60942) -
                0.0275479*(1/((0.998432 - y)**(0.492465)) - 1.90111) -
                0.157123*np.exp(-0.04738*y**8) +
                0.810239*np.exp(-40*(y - 0.133875)))
        muR = self.get_muR_from_OCV(OCV, muR_ref)
        actR = self.get_actR_None(y)
        return muR, actR

    def LiMn2O4_ss2(self, y, ybar, muR_ref, ISfuncs=None):
        """ Fuller, Doyle, Newman, 1994 """
        # OCV in V vs Li/Li+
        OCV = (4.06279 + 0.0677504*np.tanh(-21.8502*y + 12.8268) -
                0.105734*(1/((1.00167 - y)**(0.379571)) - 1.575994) -
                0.045*np.exp(-71.69*y**8) +
                0.01*np.exp(-200*(y - 0.19)))
        muR = self.get_muR_from_OCV(OCV, muR_ref)
        actR = self.get_actR_None(y)
        return muR, actR

    def LiC6_coke_ss(self, y, ybar, muR_ref, ISfuncs=None):
        """ Doyle, Newman, 1996 """
        OCV = (-0.16 + 1.32*np.exp(-3.0*y) +
                10.*np.exp(-2000.*y))
        muR = self.get_muR_from_OCV(OCV, muR_ref)
        actR = self.get_actR_None(y)
        return muR, actR

    def LiC6_coke_ss2(self, y, ybar, muR_ref, ISfuncs=None):
        """ Fuller, Doyle, Newman, 1994 """
        c1 = -0.132056
        c2 = 1.40854
        c3 = -3.52312
        OCV = c1 + c2*np.exp(c3*y)
        muR = self.get_muR_from_OCV(OCV, muR_ref)
        actR = self.get_actR_None(y)
        return muR, actR

    def LiC6_ss(self, y, ybar, muR_ref, ISfuncs=None):
        """ Safari, Delacourt 2011 """
        OCV = (0.6379 + 0.5416*np.exp(-305.5309*y) +
                0.044*np.tanh(-(y - 0.1958)/0.1088) -
                0.1978*np.tanh((y - 1.0571)/0.0854) -
                0.6875*np.tanh((y + 0.0117)/0.0529) -
                0.0175*np.tanh((y - 0.5692)/0.0875))
        muR = self.get_muR_from_OCV(OCV, muR_ref)
        actR = self.get_actR_None(y)
        return muR, actR

    def LiC6_ss2(self, y, ybar, muR_ref, ISfuncs=None):
        """ Bernardi and Go 2011 """
        p1, p2, p3, p4 = (0.085, 0.120, 0.210, 3.5)
        sfac = 0.3
        OCV = (p1*stepDown(y, 1., sfac*0.02)
                + (p2 - p1)*stepDown(y, 0.5, 0.005)
                + (p3 - p2)*stepDown(y, 0.1944, sfac*0.03571)
                + (p4 - p3)*stepDown(y, 0., sfac*0.08333))
        muR = self.get_muR_from_OCV(OCV, muR_ref)
        actR = self.get_actR_None(y)
        return muR, actR

    def LiC6_2step_ss(self, y, ybar, muR_ref, ISfuncs=None):
        """
        Fit function to the OCV predicted by the phase separating
        2-variable graphite model (LiC6 function in this class).
        """
        Vstd = 0.12
        Vstep = 0.0359646
        edgeLen = 0.024
        lEdge = edgeLen
        rEdge = 1 - edgeLen
        width = 0.001
        lSide = -((np.log(y/(1-y)) - np.log(lEdge/(1-lEdge)))
                * stepDown(y, lEdge, width))
        rSide = -((np.log(y/(1-y)) - np.log(rEdge/(1-rEdge)))
                * stepUp(y, rEdge, width))
        OCV = (
                Vstd
                + Vstep*(stepDown(y, 0.5, 0.013) - 1)
                + self.kToe*(lSide + rSide)
                )
        muR = self.get_muR_from_OCV(OCV, muR_ref)
        actR = self.get_actR_None(y)
        return muR, actR

    def Li_ss(self, y, ybar, muR_ref, ISfuncs=None):
        muR = 0.*y + muR_ref
        actR = 0.*y + 1.
        return muR, actR

    def NCA_ss1(self, y, ybar, muR_ref, ISfuncs=None):
        """
        This function was obtained from Dan Cogswell's fit of Samsung
        data.
        """
        OCV = (3.86 + 1.67*y - 9.52*y**2 + 15.04*y**3 - 7.95*y**4 -
                0.06*np.log(y/(1-y)))
        muR = self.get_muR_from_OCV(OCV, muR_ref)
        actR = self.get_actR_None(y)
        return muR, actR

    def NCA_ss2(self, y, ybar, muR_ref, ISfuncs=None):
        """
        Li_q Ni(0.8)Co(0.15)Al(0.05)O2
        as a function of y. Here, y actually represents a practical
        utilization of 70% of the material, so the material is "empty"
        (y=0) when q=0.3 and full (y=1) when q=1.
        This function was obtained from a fit by Raymond B. Smith
        of Samsung data of a LiC6-NCA cell discharged at C/100.
        """
        OCV = (-self.kToe*np.log(y/(1-y))
                + 4.12178 - 0.2338*y - 1.24566*y**2 + 1.16769*y**3
                - 0.20745*y**4)
        muR = self.get_muR_from_OCV(OCV, muR_ref)
        actR = self.get_actR_None(y)
        return muR, actR

    ######
    # Functions based on thermodynamic models
    ######

    def idealSln(self, y, ISfuncs=None):
        """ Helper function: Should not be called directly from
        simulation. Call a specific material instead. """
        T = self.T
        if ISfuncs is not None:
            # Input must be a vector when using ISfuncs
            muR = T*np.array([ISfuncs[i]() for i in range(len(y))])
        else:
            muR = T*np.log(y/(1-y))
        return muR

    def regSln(self, y, Omga, ISfuncs=None):
        """ Helper function """
        muR_IS = self.idealSln(y, ISfuncs=ISfuncs)
        enthalpyTerm = Omga*(1-2*y)
        muR = muR_IS + enthalpyTerm
        return muR

    def regSln2(self, y, Omga, Omgb, Omgc, EvdW, ISfuncs=None):
        """ Helper function """
        y1, y2 = y
        ISfuncs1, ISfuncs2 = ISfuncs
        muR1 = self.regSln(y1, Omga, ISfuncs1)
        muR2 = self.regSln(y2, Omga, ISfuncs2)
        muR1 += Omgb*y2 + Omgc*y2*(1-y2)*(1-2*y1)
        muR2 += Omgb*y1 + Omgc*y1*(1-y1)*(1-2*y2)
        muR1 += EvdW * (30 * y1**2 * (1-y1)**2)
        muR2 += EvdW * (30 * y2**2 * (1-y2)**2)
        return (muR1, muR2)

    def nonHomogRectFixedCsurf(self, y, ybar, B, kappa, ywet):
        """ Helper function """
        N = len(y)
        ytmp = np.empty(N+2, dtype=object)
        ytmp[1:-1] = y
        ytmp[0] = ywet
        ytmp[-1] = ywet
        dxs = 1./N
        curv = np.diff(ytmp, 2)/(dxs**2)
        muR_nh = -kappa*curv + B*(y - ybar)
        return muR_nh

    def nonHomogRoundWetting(self, y, ybar, B, kappa, beta_s, shape,
            r_vec):
        """ Helper function """
        dr = r_vec[1] - r_vec[0]
        Rs = 1.
        curv = calc_curv(y, dr, r_vec, Rs, beta_s, shape)
        muR_nh = B*(y - ybar) - kappa*curv
        return muR_nh

    def generalRegSln(self, y, ybar, ISfuncs):
        """ Helper function """
        ptype = self.ndD["type"]
        Omga = self.ndD["Omga"]
        N = len(y)
        muR = self.regSln(y, Omga, ISfuncs)
        if ("homog" not in ptype) and (N > 1):
            shape = self.ndD["shape"]
            kappa = self.ndD["kappa"]
            B = self.ndD["B"]
            if shape == "C3":
                cwet = self.ndD["cwet"]
                muR += self.nonHomogRectFixedCsurf(y, ybar, B, kappa, cwet)
            elif shape in ["cylinder", "sphere"]:
                beta_s = self.ndD["beta_s"]
                r_vec = mpetMaterials.get_unit_solid_discr(shape, N)[0]
                muR += self.nonHomogRoundWetting(y, ybar, B, kappa,
                        beta_s, shape, r_vec)
        return muR

    def generalRegSln2(self, y, ybar, ISfuncs):
        """ Helper function """
        ptype = self.ndD["type"]
        Omga = self.ndD["Omga"]
        Omgb = self.ndD["Omgb"]
        Omgc = self.ndD["Omgc"]
        EvdW = self.ndD["EvdW"]
        N = len(y[0])
        muR1, muR2 = self.regSln2(y, Omga, Omgb, Omgc, EvdW, ISfuncs)
        if ("homog" not in ptype) and (N > 1):
            shape = self.ndD["shape"]
            kappa = self.ndD["kappa"]
            B = self.ndD["B"]
            beta_s = self.ndD["beta_s"]
            r_vec = mpetMaterials.get_unit_solid_discr(shape, N)[0]
            muR1 += self.nonHomogRoundWetting(y[0], ybar[0], B,
                    kappa, beta_s, shape, r_vec)
            muR2 += self.nonHomogRoundWetting(y[1], ybar[1], B,
                    kappa, beta_s, shape, r_vec)
        return (muR1, muR2)

    def LiFePO4(self, y, ybar, muR_ref, ISfuncs=None):
        """ Bai, Cogswell, Bazant 2011 """
        muRtheta = -self.eokT*3.422
        muR = self.generalRegSln(y, ybar, ISfuncs)
        actR = np.exp(muR/self.T)
        muR += muRtheta + muR_ref
        return muR, actR

    def LiC6(self, y, ybar, muR_ref, ISfuncs=(None, None)):
        """ Ferguson and Bazant 2014 """
        muRtheta = -self.eokT*0.12
        muR1, muR2 = self.generalRegSln2(y, ybar, ISfuncs)
        actR1 = np.exp(muR1/self.T)
        actR2 = np.exp(muR2/self.T)
        muR1 += muRtheta + muR_ref
        muR2 += muRtheta + muR_ref
        return (muR1, muR2), (actR1, actR2)

def stepDown(x, xc, delta):
    return 0.5*(-np.tanh((x - xc)/delta) + 1)
def stepUp(x, xc, delta):
    return 0.5*(np.tanh((x - xc)/delta) + 1)

def calc_curv(c, dr, r_vec, Rs, beta_s, particleShape):
    N = len(c)
    curv = np.empty(N, dtype=object)
    if particleShape == "sphere":
        curv[0] = 3 * (2*c[1] - 2*c[0]) / dr**2
        curv[1:N-1] = (np.diff(c, 2)/dr**2 +
                (c[2:] - c[0:-2])/(dr*r_vec[1:-1]))
        curv[N-1] = ((2./Rs)*beta_s +
                (2*c[-2] - 2*c[-1] + 2*dr*beta_s)/dr**2)
    elif particleShape == "cylinder":
        curv[0] = 2 * (2*c[1] - 2*c[0]) / dr**2
        curv[1:N-1] = (np.diff(c, 2)/dr**2 +
                (c[2:] - c[0:-2])/(2 * dr*r_vec[1:-1]))
        curv[N-1] = ((1./Rs)*beta_s +
                (2*c[-2] - 2*c[-1] + 2*dr*beta_s)/dr**2)
    else:
        raise NotImplementedError("calc_curv_c only for sphere and cylinder")
    return curv

