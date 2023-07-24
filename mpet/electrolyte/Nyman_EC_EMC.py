from mpet.config import constants

def Nyman_EC_EMC():
    """ Set of parameters from Nyman, Behm, Lindberg 2007,
    for LiPF6 in EC-EMC 3:7
    """

    def tp0(c, T):
        ce = c*1  # dimensionalized c
        (p0, p1, p2, 
         p3, ) = (
            0.4492, -0.4717, 0.4106, 
            -0.1287)
        return p0 + p1*ce + p2*ce**2 + p3*ce**3 

    def D(c, T):
        ce = c*1  # dimensionalized c
        (p0, p1, p2) = (4.862*1e-10, -3.972*1e-10, 8.794*1e-11)
        return p0 + p1*ce + p2*ce**2 # m^2/s
    
    def therm_fac(c, T):
        ce = c*1  # dimensionalized c
        (p0, p1, p2) = (0.44103, -0.74678, 0.28687)
        tmp = p0 + p1*ce + p2*ce**2
        return tmp/(1-tp0(c, T))

    def sigma(c, T):
        ce = c*1  # dimensionalized c
        (A1, A2, A3) = (1.297, -25.1, 33.29)
        out = A1*ce**3 + A2*ce**1.5 + A3*ce # mS/cm
        out = out*0.1 # S/m   
        return out
    
    Dref = D(constants.c_ref/1000, 1)

    def D_ndim(c, T):
        return D(c, T) / Dref

    def sigma_ndim(c, T):
        return sigma(c, T) * (
            constants.k*constants.T_ref/(constants.e**2*Dref*constants.N_A*(constants.c_ref)))
    
    return D_ndim, sigma_ndim, therm_fac, tp0, Dref