"""
Manage functions for the parameters involved in Stefan-Maxwell
based concentrated electrolyte transport theory for binary
electrolytes.

Each electrolyte set must output functions for the following as a
function of c (electrolyte concentration, M)
 - Dchem [m^2/s] = the prefactor for grad(c) in species conservation
 - kappa [S/m] = the conductivity
 - (1 + dln(f_\pm)/dln(c)) = the "thermodynamic factor"
 - t_+^0 = the transference number of the cations
"""

cref = 1. # M
Tref = 298. # K
N_A = 6.022e23
k = 1.381e-23
e = 1.602e-19

def test1():
#    cref = self.cref
#    Tref = self.Tref
#    N_A = self.N_A
#    k = self.k
#    e = self.e
    """ Set of dilute solution parameters with zp=|zm|=nup=num=1,
    Dp = 2.2e-10 m^2/s
    Dm = 2.94e-10 m^2/s
    """
    Dp = 2.2e-10
    Dm = 2.94e-10
    def D(c):
        return (2*Dp*Dm/(Dp+Dm))
    def thermFac(c):
        return 1.
    def tp0(c):
        return Dp/(Dp+Dm)
    def kappa(c):
        return Dm*c*N_A*e**2/(k*Tref*(1-tp0(c)))
    Dref = D(cref)
    D_ndim = lambda c: D(c) / Dref
    kappa_ndim = lambda c: kappa(c) * (
            k*Tref/(e**2*Dref*N_A*cref))
    return D_ndim, kappa_ndim, thermFac, tp0, Dref

def getProps(elytePropSet):
    if elytePropSet == "test1":
        return test1()
    else:
        raise NotImplementedError("Unrecognized elytePropSet")
