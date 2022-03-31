"""This module handles properties associated with the active materials."""
import numpy as np

import mpet.geometry as geo
from mpet.config import constants


class Dfuncs():
    """This class returns the filling-fraction dependent variation of
    the transport coefficient, D(y), such that
    Flux = -D_ref*D(y)*grad(y) for solid solution transport or
    Flux = -D_ref*D(y)*grad(mu) for thermo-based transport
    where y here is the filling fraction, D_ref has dimensions of
    length^2/time, D(y) is dimensionless, and mu, the chemical
    potential, has been scaled to the thermal energy, k*Tref. For more
    details on the two forms, see muRfuncs which defines both solid
    solution (_ss) and materials based on simpler thermodynamic models.
    """

    def __init__(self, Dfunc):
        Dopts = {}
        Dopts['lattice'] = self.lattice
        Dopts['constant'] = self.constant
        Dopts['NMC532_Colclasure20'] = self.NMC532_Colclasure20
        self.Dfunc = Dopts[Dfunc]

    def constant(self, y):
        return 1.

    def lattice(self, y):
        return y*(1-y)

    def NMC532_Colclasure20(self,y):
        cm2_s = 20000.*10**(-2319.*y**10 + 6642.*y**9 - 5269.*y**8 - 3319.*y**7 + 10038.*y**6
                            - 9806.*y**5 + 5817.*y**4 - 2286.*y**3 + 575.3*y**2 - 83.16*y-9.292)
        return(cm2_s/1e4)


class muRfuncs():
    """ This class defines functions which describe the chemical
    potential of active materials.
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
    def __init__(self, config, trode, ind=None):
        """config is the full dictionary of
        parameters for the electrode particles, as made for the
        simulations. trode is the selected electrode.
        ind is optinally the selected particle, provided as (vInd, pInd)
        """
        self.config = config
        self.trode = trode
        self.ind = ind
        self.T = config['T']  # nondimensional
        # eokT and kToe are the reference values for scalings
        self.eokT = constants.e / (constants.k * constants.T_ref)
        self.kToe = 1. / self.eokT

        # Convert "muRfunc" to a callable function
        self.muRfunc = getattr(self, self.get_trode_param("muRfunc"))

    def get_trode_param(self, item):
        """
        Shorthand to retrieve electrode-specific value
        """
        value = self.config[self.trode, item]
        # check if it is a particle-specific parameter
        if self.ind is not None and item in self.config.params_per_particle:
            value = value[self.ind]
        return value

    def get_muR_from_OCV(self, OCV, muR_ref):
        return -self.eokT*OCV + muR_ref

    ######
    # Solid solution functions
    # These are all obtained directly from fitting an OCV of the
    # material with a standard counter electrode.
    # They can all only return values at 298 K
    ######

    def LiMn2O4_ss(self, y, ybar, muR_ref, ISfuncs=None):
        """ Doyle, Newman, 1996 """
        # OCV in V vs Li/Li+
        OCV = (4.19829 + 0.0565661*np.tanh(-14.5546*y + 8.60942)
               - 0.0275479*(1/((0.998432 - y)**(0.492465)) - 1.90111)
               - 0.157123*np.exp(-0.04738*y**8)
               + 0.810239*np.exp(-40*(y - 0.133875)))
        muR = self.get_muR_from_OCV(OCV, muR_ref)
        actR = None
        return muR, actR

    def LiMn2O4_ss2(self, y, ybar, muR_ref, ISfuncs=None):
        """ Fuller, Doyle, Newman, 1994 """
        # OCV in V vs Li/Li+
        OCV = (4.06279 + 0.0677504*np.tanh(-21.8502*y + 12.8268)
               - 0.105734*(1/((1.00167 - y)**(0.379571)) - 1.575994)
               - 0.045*np.exp(-71.69*y**8)
               + 0.01*np.exp(-200*(y - 0.19)))
        muR = self.get_muR_from_OCV(OCV, muR_ref)
        actR = None
        return muR, actR

    def LiC6_coke_ss(self, y, ybar, muR_ref, ISfuncs=None):
        """ Doyle, Newman, 1996 """
        OCV = (-0.16 + 1.32*np.exp(-3.0*y) + 10.*np.exp(-2000.*y))
        muR = self.get_muR_from_OCV(OCV, muR_ref)
        actR = None
        return muR, actR

    def LiC6_coke_ss2(self, y, ybar, muR_ref, ISfuncs=None):
        """ Fuller, Doyle, Newman, 1994 """
        c1 = -0.132056
        c2 = 1.40854
        c3 = -3.52312
        OCV = c1 + c2*np.exp(c3*y)
        muR = self.get_muR_from_OCV(OCV, muR_ref)
        actR = None
        return muR, actR

    def LiC6_ss(self, y, ybar, muR_ref, ISfuncs=None):
        """ Safari, Delacourt 2011 """
        OCV = (0.6379 + 0.5416*np.exp(-305.5309*y)
               + 0.044*np.tanh(-(y - 0.1958)/0.1088)
               - 0.1978*np.tanh((y - 1.0571)/0.0854)
               - 0.6875*np.tanh((y + 0.0117)/0.0529)
               - 0.0175*np.tanh((y - 0.5692)/0.0875))
        muR = self.get_muR_from_OCV(OCV, muR_ref)
        actR = None
        return muR, actR

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

    def LiC6_2step_ss(self, y, ybar, muR_ref, ISfuncs=None):
        """
        Fit function to the OCV predicted by the phase separating
        2-variable graphite model (LiC6 function in this class).
        """
        Vstd = 0.1196
        Vstep = 0.0351733976
        edgeLen = 0.024
        lEdge = edgeLen
        rEdge = 1 - edgeLen
        width = 1e-4
        vshift = 1e-2
        lSide = -((np.log(y/(1-y)) - np.log(lEdge/(1-lEdge)) - vshift)*step_down(y, lEdge, width))
        rSide = -((np.log(y/(1-y)) - np.log(rEdge/(1-rEdge)) + vshift)*step_up(y, rEdge, width))
        OCV = (
            Vstd
            + Vstep*(step_down(y, 0.5, 0.013) - 1)
            + self.kToe*(lSide + rSide)
            )
        muR = self.get_muR_from_OCV(OCV, muR_ref)
        actR = None
        return muR, actR

    def LiC6_Colclasure_1506T(self, y, ybar, muR_ref, ISfuncs=None):
        """Taken from Colclasure 2020, but only for lithiating. The range of\
        confidence is 0.01 to 0.97 for voltages of ~0.6V and \
        0.0435V. Don’t calculate U2 for intercalation fractions below 0.8, \
        it will cause numeric issues just use U1."""
        x = y
        U1 = \
            - 1.059423355572770E-02*np.tanh((x - 1.453708425609560E-02) / 9.089868397988610E-05) \
            + 2.443615203087110E-02*np.tanh((x - 5.464261369950400E-01) / 6.270508166379020E-01) \
            - 1.637520788053810E-02*np.tanh((x - 5.639025014475490E-01) / 7.053886409518520E-02) \
            - 6.542365622896410E-02*np.tanh((x - 5.960370524233590E-01) / 1.409966536648620E+00) \
            - 4.173226059293490E-02*np.tanh((x - 1.787670587868640E-01) / 7.693844911793470E-02) \
            - 4.792178163846890E-01*np.tanh((x + 3.845707852011820E-03) / 4.112633446959460E-02) \
            + 6.594735004847470E-01 \
            - 4.364293924074990E-02*np.tanh((x - 9.449231893318330E-02) / -2.046776012570780E-02) \
            - 8.241166396760410E-02*np.tanh((x - 7.746685789572230E-02) / 3.593817905677970E-02)
        U2 = \
            -1.731504647676420E+02*np.power(x, 8.0) + 8.252008712749000E+01*np.power(x, 7.0) \
            + 1.233160814852810E+02*np.power(x, 6.0) + 5.913206621637760E+01*np.power(x, 5.0) \
            + 3.322960033709470E+01*np.power(x, 4.0) + 3.437968012320620E+00*np.power(x, 3.0) \
            - 6.906367679257650E+01*np.power(x, 2.0) - \
            1.228217254296760E+01*x - 5.037944982759270E+01
        U = U1 + (U2 - U1) / (1.0 + np.exp(-1.0e2*(x - 1.02956203215198)))
        muR = self.get_muR_from_OCV(U, muR_ref)
        actR = None
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
        OCV = (3.86 + 1.67*y - 9.52*y**2 + 15.04*y**3 - 7.95*y**4
               - 0.06*np.log(y/(1-y)))
        muR = self.get_muR_from_OCV(OCV, muR_ref)
        actR = None
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
        actR = None
        return muR, actR

    def NMC532_Colclasure20(self, y, ybar, muR_ref, ISfuncs=None):
        x = y
        OCV = (5.314735633000300E+00 +
               -3.640117692001490E+03*x**14.0 + 1.317657544484270E+04*x**13.0
               - 1.455742062291360E+04*x**12.0 - 1.571094264365090E+03*x**11.0
               + 1.265630978512400E+04*x**10.0 - 2.057808873526350E+03*x**9.0
               - 1.074374333186190E+04*x**8.0 + 8.698112755348720E+03*x**7.0
               - 8.297904604107030E+02*x**6.0 - 2.073765547574810E+03*x**5.0
               + 1.190223421193310E+03*x**4.0 - 2.724851668445780E+02*x**3.0
               + 2.723409218042130E+01*x**2.0 - 4.158276603609060E+00*x +
               -5.573191762723310E-04*np.exp(6.560240842659690E+00*x**4.148209275061330E+01)
               )
        muR = self.get_muR_from_OCV(OCV, muR_ref)
        actR = None
        return muR, actR

    def LTO(self, y, ybar, muR_ref, ISfuncs=None):
        """
        Vasileiadis 2017
        """
        muRtheta = -self.eokT * 1.55
        muRhomog = self.reg_sln(y, self.get_trode_param("Omega_a"), ISfuncs)
        muRnonHomog = self.general_non_homog(y, ybar)
        muR = muRhomog + muRnonHomog
        actR = np.exp(muR / self.T)
        muR += muRtheta + muR_ref
        return muR, actR

    def testIS_ss(self, y, ybar, muR_ref, ISfuncs=None):
        """Ideal solution material for testing."""
        OCV = -self.kToe*np.log(y/(1-y))
        muR = self.get_muR_from_OCV(OCV, muR_ref)
        actR = None
        return muR, actR

    def testRS_ss(self, y, ybar, muR_ref, ISfuncs=None):
        """
        Regular solution material which phase separates at binodal points,
        for modeling as a solid solution. For testing.
        """
        # Based Omg = 3*k*T_ref
        yL = 0.07072018
        yR = 0.92927982
        OCV_rs = -self.kToe*self.reg_sln(y, self.get_trode_param("Omega_a"), ISfuncs)
        width = 0.005
        OCV = OCV_rs*step_down(y, yL, width) + OCV_rs*step_up(y, yR, width) + 2
        muR = self.get_muR_from_OCV(OCV, muR_ref)
        actR = None
        return muR, actR

    ######
    # Functions based on thermodynamic models
    ######

    def ideal_sln(self, y, ISfuncs=None):
        """ Helper function: Should not be called directly from
        simulation. Call a specific material instead. """
        T = self.T
        if ISfuncs is not None:
            # Input must be a vector when using ISfuncs
            muR = T*np.array([ISfuncs[i]() for i in range(len(y))])
        else:
            muR = T*np.log(y/(1-y))
        return muR

    def reg_sln(self, y, Omga, ISfuncs=None):
        """ Helper function """
        muR_IS = self.ideal_sln(y, ISfuncs=ISfuncs)
        enthalpyTerm = Omga*(1-2*y)
        muR = muR_IS + enthalpyTerm
        return muR

    def graphite_2param_homog(self, y, Omga, Omgb, Omgc, EvdW, ISfuncs=None):
        """ Helper function """
        y1, y2 = y
        if ISfuncs is None:
            ISfuncs1, ISfuncs2 = None, None
        else:
            ISfuncs1, ISfuncs2 = ISfuncs
        muR1 = self.reg_sln(y1, Omga, ISfuncs1)
        muR2 = self.reg_sln(y2, Omga, ISfuncs2)
        muR1 += Omgb*y2 + Omgc*y2*(1-y2)*(1-2*y1)
        muR2 += Omgb*y1 + Omgc*y1*(1-y1)*(1-2*y2)
        muR1 += EvdW * (30 * y1**2 * (1-y1)**2)
        muR2 += EvdW * (30 * y2**2 * (1-y2)**2)
        return (muR1, muR2)

    def graphite_1param_homog(self, y, Omga, Omgb, ISfuncs=None):
        """ Helper function """
        width = 5e-2
        tailScl = 5e-2
        muLtail = -tailScl*1./(y**(0.85))
        muRtail = tailScl*1./((1-y)**(0.85))
        slpScl = 0.45
        muLlin = slpScl*Omga*4*(0.26-y)*step_down(y, 0.5, width)
        muRlin = (slpScl*Omga*4*(0.74-y) + Omgb)*step_up(y, 0.5, width)
        muR = muLtail + muRtail + muLlin + muRlin
        return muR

    def graphite_1param_homog_2(self, y, Omga, Omgb, ISfuncs=None):
        """ Helper function """
        width = 5e-2
        tailScl = 5e-2
        slpScl = 0.45
        muLtail = -tailScl*1./(y**(0.85))
        muRtail = tailScl*1./((1-y)**(0.85))
        muLlin = (slpScl*Omga*12*(0.40-y)
                  * step_down(y, 0.49, 0.9*width)*step_up(y, 0.35, width))
        muRlin = (slpScl*Omga*4*(0.74-y) + Omgb)*step_up(y, 0.5, width)
        muLMod = (0.
                  + 40*(-np.exp(-y/0.015))
                  + 0.75*(np.tanh((y-0.17)/0.02) - 1)
                  + 1.0*(np.tanh((y-0.22)/0.040) - 1)
                  )*step_down(y, 0.35, width)
        muR = muLMod + muLtail + muRtail + muLlin + muRlin
        return muR

    def graphite_1param_homog_3(self, y, Omga, Omgb, ISfuncs=None):
        """ Helper function with low hysteresis and soft tail """
        width = 5e-2
        tailScl = 5e-2
        muLtail = -tailScl*1./(y**(0.85))
        muRtail = tailScl*1./((1-y)**(0.85))
        muRtail = 1.0e1*step_up(y, 1.0, 0.045)
        muLlin = (0.15*Omga*12*(0.40-y**0.98)
                  * step_down(y, 0.49, 0.9*width)*step_up(y, 0.35, width))
        muRlin = (0.1*Omga*4*(0.74-y) + 0.90*Omgb)*step_up(y, 0.5, 0.4*width)
        muLMod = (0.
                  + 40*(-np.exp(-y/0.015))
                  + 0.75*(np.tanh((y-0.17)/0.02) - 1)
                  + 1.0*(np.tanh((y-0.22)/0.040) - 1)
                  )*step_down(y, 0.35, width)
        muR = 0.18 + muLMod + muLtail + muRtail + muLlin + muRlin
        return muR

    def non_homog_rect_fixed_csurf(self, y, ybar, B, kappa, ywet):
        """ Helper function """
        if isinstance(y, np.ndarray):
            N = len(y)
            ytmp = np.empty(N+2, dtype=object)
            ytmp[1:-1] = y
            ytmp[0] = ywet
            ytmp[-1] = ywet
            dxs = 1./N
            curv = np.diff(ytmp, 2)/(dxs**2)
            muR_nh = -kappa*curv + B*(y - ybar)
        elif (isinstance(y, tuple) and len(y) == 2
                and isinstance(y[0], np.ndarray)):
            stech_Mn = 0.5
            stech_Fe = 0.5
            ybar_avg = stech_Mn*ybar[0]+stech_Fe*ybar[1]
            N = len(y[0])

            ytmp1 = np.empty(N+2, dtype=object)
            ytmp1[1:-1] = y[0]
            ytmp1[0] = ywet
            ytmp1[-1] = ywet
            dxs = 1./N
            curv1 = np.diff(ytmp1, 2)/(dxs**2)
            muR1_nh = (-kappa*curv1
                       + B*(y[0] - ybar[0]))

            ytmp2 = np.empty(N+2, dtype=object)
            ytmp2[1:-1] = y[1]
            ytmp2[0] = ywet
            ytmp2[-1] = ywet
            dxs = 1./N
            curv2 = np.diff(ytmp2, 2)/(dxs**2)
            muR2_nh = (-kappa*curv2 
                       + B*(y[1] - ybar[1]))

            muR_nh = (muR1_nh, muR2_nh)
        return muR_nh
    
    def non_homog_rect_variational(self, y, ybar, B, kappa):
        """ Helper function """
        # the taylor expansion at the edges is used
        if isinstance(y, np.ndarray):
            N_2 = len(y)
            ytmp = np.empty(N_2+2, dtype=object)
            dxs = 1./N_2
            ytmp[1:-1] = y
            ytmp[0] = y[0] + np.diff(y)[0]*dxs + 0.5*np.diff(y,2)[0]*dxs**2
            ytmp[-1] = y[-1] + np.diff(y)[-1]*dxs + 0.5*np.diff(y,2)[-1]*dxs**2
            curv = np.diff(ytmp, 2)/(dxs**2)
            muR_nh = -kappa*curv + B*(y - ybar)
        elif (isinstance(y, tuple) and len(y) == 2
                and isinstance(y[0], np.ndarray)):
            y1 = y[0]
            y2 = y[1]
            stech_Mn = 0.5
            stech_Fe = 0.5
            ybar_avg = stech_Mn*ybar[0]+stech_Fe*ybar[1]
            N_2 = len(y[0])

            ytmp1 = np.empty(N_2+2, dtype=object)
            dxs = 1./N_2
            ytmp1[1:-1] = y1
            ytmp1[0] = y1[0] + np.diff(y1)[0]*dxs + 0.5*np.diff(y1,2)[0]*dxs**2
            ytmp1[-1] = y1[-1] + np.diff(y1)[-1]*dxs + 0.5*np.diff(y1,2)[-1]*dxs**2
            curv1 = np.diff(ytmp1, 2)/(dxs**2)
            muR1_nh = -kappa*curv1 + B*(y1 - ybar[0])

            ytmp2 = np.empty(N_2+2, dtype=object)
            dxs = 1./N_2
            ytmp2[1:-1] = y2
            ytmp2[0] = y2[0] + np.diff(y2)[0]*dxs + 0.5*np.diff(y2,2)[0]*dxs**2
            ytmp2[-1] = y2[-1] + np.diff(y2)[-1]*dxs + 0.5*np.diff(y2,2)[-1]*dxs**2
            curv2 = np.diff(ytmp2, 2)/(dxs**2)
            muR2_nh = -kappa*curv2 + B*(y2 - ybar[1])

            muR_nh = (muR1_nh, muR2_nh)
        return muR_nh

    def non_homog_round_wetting(self, y, ybar, B, kappa, beta_s, shape, r_vec):
        """ Helper function """
        dr = r_vec[1] - r_vec[0]
        Rs = 1.
        curv = geo.calc_curv(y, dr, r_vec, Rs, beta_s, shape)
        muR_nh = B*(y - ybar) - kappa*curv
        return muR_nh

    def general_non_homog(self, y, ybar):
        """ Helper function """
        ptype = self.get_trode_param("type")
        mod1var, mod2var = False, False
        if isinstance(y, np.ndarray):
            mod1var = True
            N = len(y)
        elif (isinstance(y, tuple) and len(y) == 2
                and isinstance(y[0], np.ndarray)):
            mod2var = True
            N = len(y[0])
        else:
            raise Exception("Unknown input type")
        if ("homog" not in ptype) and (N > 1):
            shape = self.get_trode_param("shape")
            if shape == "C3":
                if mod1var:
                    kappa = self.get_trode_param("kappa")
                    B = self.get_trode_param("B")
                    if self.get_trode_param("type") in ["ACR"]:
                        cwet = self.get_trode_param("cwet")
                        muR_nh = self.non_homog_rect_fixed_csurf(
                            y, ybar, B, kappa, cwet)
                    elif self.get_trode_param("type") in ["ACR_Diff"]:
                        muR_nh = self.non_homog_rect_variational(
                            y, ybar, B, kappa)
                elif mod2var:
                    kappa= self.get_trode_param("kappa")
                    B = self.get_trode_param("B")
                    cwet = self.get_trode_param("cwet")
                    muR_nh = self.non_homog_rect_fixed_csurf(
                        y, ybar, B, kappa, cwet)
                    # muR_nh = self.non_homog_rect_variational(
                    #     y, ybar, B, kappa)
            elif shape in ["cylinder", "sphere"]:
                kappa = self.get_trode_param("kappa")
                B = self.get_trode_param("B")
                beta_s = self.get_trode_param("beta_s")
                r_vec = geo.get_unit_solid_discr(shape, N)[0]
                if mod1var:
                    muR_nh = self.non_homog_round_wetting(
                        y, ybar, B, kappa, beta_s, shape, r_vec)
                elif mod2var:
                    muR1_nh = self.non_homog_round_wetting(
                        y[0], ybar[0], B, kappa, beta_s, shape, r_vec)
                    muR2_nh = self.non_homog_round_wetting(
                        y[1], ybar[1], B, kappa, beta_s, shape, r_vec)
                    muR_nh = (muR1_nh, muR2_nh)
        else:  # homogeneous particle
            if mod1var:
                muR_nh = 0*y
            elif mod2var:
                muR_nh = (0*y[0], 0*y[1])
        return muR_nh

    def LiFePO4(self, y, ybar, muR_ref, ISfuncs=None):
        """ Bai, Cogswell, Bazant 2011 """
        muRtheta = -self.eokT*3.422
        muRhomog = self.reg_sln(y, self.get_trode_param("Omega_a"), ISfuncs)
        muRnonHomog = self.general_non_homog(y, ybar)
        muR = muRhomog + muRnonHomog
        actR = np.exp(muR/self.T)
        muR += muRtheta + muR_ref
        return muR, actR
        
    def LiFeMnPO4(self, y, ybar, muR_ref, ISfuncs = None):
        """ New test material """
        muRtheta1 = -self.eokT*4.1
        muRtheta2 = -self.eokT*3.9
        if ISfuncs is None:
            ISfuncs1, ISfuncs2 = None, None
        else:
            ISfuncs1, ISfuncs2 = ISfuncs
        y1, y2 = y
        Omgc = self.get_trode_param('Omega_c')
        muR1 = self.reg_sln(y1, self.get_trode_param('Omega_a'), ISfuncs1)
        muR2 = self.reg_sln(y2, self.get_trode_param('Omega_b'), ISfuncs2)
        muR1 += Omgc*y2*y1
        muR2 += Omgc*y1*y2
        muR1nonHomog, muR2nonHomog = self.general_non_homog(y, ybar)
        muR1 += muR1nonHomog
        muR2 += muR2nonHomog
        actR1 = np.exp(muR1/self.T)
        actR2 = np.exp(muR2/self.T)
        muR1 += muRtheta1 + muR_ref
        muR2 += muRtheta2 + muR_ref
        return (muR1, muR2), (actR1, actR2)

    def LiC6(self, y, ybar, muR_ref, ISfuncs=(None, None)):
        """ Ferguson and Bazant 2014 """
        muRtheta = -self.eokT*0.12
        muR1homog, muR2homog = self.graphite_2param_homog(
            y, self.get_trode_param("Omega_a"), self.get_trode_param("Omega_b"),
            self.get_trode_param("Omega_c"), self.get_trode_param("EvdW"), ISfuncs)
        muR1nonHomog, muR2nonHomog = self.general_non_homog(y, ybar)
        muR1 = muR1homog + muR1nonHomog
        muR2 = muR2homog + muR2nonHomog
        actR1 = np.exp(muR1/self.T)
        actR2 = np.exp(muR2/self.T)
        muR1 += muRtheta + muR_ref
        muR2 += muRtheta + muR_ref
        return (muR1, muR2), (actR1, actR2)

    def LiC6_1param(self, y, ybar, muR_ref, ISfuncs=None):
        muRtheta = -self.eokT*0.12
        muRhomog = self.graphite_1param_homog_3(
            y, self.get_trode_param("Omega_a"), self.get_trode_param("Omega_b"), ISfuncs)
        muRnonHomog = self.general_non_homog(y, ybar)
        muR = muRhomog + muRnonHomog
        actR = np.exp(muR/self.T)
        muR += muRtheta + muR_ref
        return muR, actR

    def testRS(self, y, ybar, muR_ref, ISfuncs=None):
        muRtheta = 0.
        muR = self.reg_sln(y, self.get_trode_param("Omega_a"), ISfuncs)
        actR = np.exp(muR/self.T)
        muR += muRtheta + muR_ref
        return muR, actR

    def testRS_ps(self, y, ybar, muR_ref, ISfuncs=None):
        muRtheta = -self.eokT*2.
        muRhomog = self.reg_sln(y, self.get_trode_param("Omega_a"), ISfuncs)
        muRnonHomog = self.general_non_homog(y, ybar)
        muR = muRhomog + muRnonHomog
        actR = np.exp(muR/self.T)
        muR += muRtheta + muR_ref
        return muR, actR

    def LiCoO2_LIONSIMBA(self, y, ybar, muR_ref, ISfuncs=None):
        """ Torchio et al, 2016. """
        T = self.T
        r1 = 4.656
        r2 = 88.669
        r3 = 401.119
        r4 = 342.909
        r5 = 462.471
        r6 = 433.434
        r7 = 1
        r8 = 18.933
        r9 = 79.532
        r10 = 37.311
        r11 = 73.083
        r12 = 95.96
        OCV_ref = (-r1 + r2*y**2 - r3*y**4 + r4*y**6 - r5*y**8 + r6 * y**10) / \
            (-r7 + r8*y**2 - r9*y**4 + r10*y**6 - r11*y**8 + r12*y**10)
        k1 = -0.001
        k2 = 0.199521039
        k3 = -0.928373822
        k4 = 1.364550689000003
        k5 = -0.6115448939999998
        k6 = 1
        k7 = -5.661479886999997
        k8 = 11.47636191
        k9 = -9.82431213599998
        k10 = 3.048755063
        dUdT = k1*(k2+k3*y+k4*y**2+k5*y**3)/(k6+k7*y+k8*y**2+k9*y**3+k10*y**4)
        OCV = OCV_ref + dUdT*(T-1)*constants.T_ref
        muR = self.get_muR_from_OCV(OCV, muR_ref)
        actR = None
        return muR, actR

    def LiC6_LIONSIMBA(self, y, ybar, muR_ref, ISfuncs=None):
        """ Torchio et al, 2016. """
        T = self.T
        r1 = 0.7222
        r2 = 0.1387
        r3 = 0.029
        r4 = 0.0172
        r5 = 0.0019
        r6 = 0.2808
        r7 = 0.7984
        OCV_ref = r1 + r2*y + r3*y**0.5 - r4 * \
            y**(-1) + r5*y**(-1.5) + r6*np.exp(0.9-15*y) - r7*np.exp(0.4465*y-0.4108)
        k1 = 0.001
        k2 = 0.005269056
        k3 = 3.299265709
        k4 = -91.79325798
        k5 = 1004.911008
        k6 = -5812.278127
        k7 = 19329.7549
        k8 = -37147.8947
        k9 = 38379.18127
        k10 = -16515.05308
        k11 = 1
        k12 = -48.09287227
        k13 = 1017.234804
        k14 = -10481.80419
        k15 = 59431.3
        k16 = -195881.6488
        k17 = 374577.3152
        k18 = -385821.1607
        k19 = 165705.8597
        dUdT = k1*(k2+k3*y+k4*y**2+k5*y**3+k6*y**4+k7*y**5+k8*y**6+k9*y**7+k10*y**8) / \
            (k11+k12*y+k13*y**2+k14*y**3+k15*y**4+k16*y**5+k17*y**6+k18*y**7+k19*y**8)
        OCV = OCV_ref + dUdT*(T-1)*constants.T_ref
        muR = self.get_muR_from_OCV(OCV, muR_ref)
        actR = None
        return muR, actR


def step_down(x, xc, delta):
    return 0.5*(-np.tanh((x - xc)/delta) + 1)


def step_up(x, xc, delta):
    return 0.5*(np.tanh((x - xc)/delta) + 1)
