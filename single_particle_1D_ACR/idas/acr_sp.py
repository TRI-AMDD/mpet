#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import sys
import os
from time import localtime, strftime
import math

import numpy as np
import scipy.sparse as sprs
import scipy.interpolate as sint

from daetools.pyDAE import *
from daetools.pyDAE.data_reporters import *
from daetools.solvers.superlu import pySuperLU
#from daetools.solvers.superlu_mt import pySuperLU_MT
from daetools.solvers.trilinos import pyTrilinos
#from daetools.solvers.intel_pardiso import pyIntelPardiso
from pyUnits import s
#from pyUnits import s, kg, m, K, Pa, mol, J, W

#########################################################################
# CONSTANTS
k = 1.381e-23      # Boltzmann constant
T = 298            # Temp, K
Tref = 298         # Reference temp, K (for non-dimensionalization)
e = 1.602e-19      # Charge of proton, C
Na = 6.02e23       # Avogadro's number
F = e*Na           # Faraday's number

#########################################################################
# SET DIMENSIONAL VALUES HERE
#init_voltage = 3.5                 # Initial cell voltage, V
currset = 0.001                         # Battery discharge c-rate
partsize = 50                       # particle size, nm

# Particle size
part_size = partsize * 1e-9                  # Average particle size, m

# Material properties
dim_a = 1.8560e-20                 # Regular solution parameter, J
dim_kappa = 5.0148e-10             # Gradient penalty, J/m
dim_b = 0.1916e9                   # Stress, Pa
rhos = 1.3793e28                   # site density, 1/m^3
csmax = rhos/Na                    # maximum concentration, mol/m^3
cwet = 0.98                        # Dimensionless wetted conc.
wet_thick = 2e-9                   # Thickness of wetting on surf.
Vstd = 3.422                       # Standard potential, V
alpha = 0.5                        # Charge transfer coefficient

# Discretization settings
solid_disc = 1e-9                 # Discretization size of solid, m (MUST BE LESS THAN LAMBDA)
tsteps = 200                      # Number disc. in time


# Define some variable types
mole_frac_t = daeVariableType(name="mole_frac_t", units=unit(),
#        lowerBound=-1e+20, upperBound=1e+20, initialGuess=0,
        lowerBound=0, upperBound=1, initialGuess=0.25,
        absTolerance=1e-6)
elec_pot_t = daeVariableType(name="elec_pot_t", units=unit(),
        lowerBound=-1e20, upperBound=1e20, initialGuess=0,
        absTolerance=1e-5)

class noise(daeScalarExternalFunction):
    def __init__(self, Name, Model, units, time, noise_vec, numnoise,
            currset):
        arguments = {}
        tmax = 1/currset.GetValue()
        self.counter = 0
        self.previous_value = None
        self.interp = sint.interp1d(
                tmax*np.arange(numnoise), noise_vec, axis=0)
        arguments["time"] = time
        daeScalarExternalFunction.__init__(self, Name, Model, units, arguments)

    def Calculate(self, values):
        time = values["time"]
        # A derivative for Jacobian is requested - return always 0.0
        if time.Derivative != 0:
            return adouble(0)
        # Store the previous time value to prevent excessive
        # interpolation.
        if self.previous_value and self.previous_value[0] == time.Value:
            return self.previous_value[1]
        # interp returns ndarrays, in this case 0-dimensional
        # therefore, the values should be converted to float
        noise_val = float(self.interp(time.Value))
        self.counter += 1
        self.previous_value = (time.Value, adouble(noise_val))
        return adouble(noise_val)

class modRay(daeModel):
    def __init__(self, Name, Parent = None, Description = ""):
        print "Init modRay"
        daeModel.__init__(self, Name, Parent, Description)

        # Domains where variables are distributed
        self.N = daeDomain("N", self, unit(), "")
        
        # variables distributed on N
        self.c = daeVariable("c", mole_frac_t, self,
                "Concentration", [self.N])
        self.cbar = daeVariable("cbar", mole_frac_t, self,
                "Average concentration in the particle")
        self.phi = daeVariable("phi", elec_pot_t, self,
                "Electrostatic potential in the solid")

        # Parameters
        self.a = daeParameter("a", unit(), self,
                "regular solution parameter")
        self.b = daeParameter("b", unit(), self,
                "Stress")
        self.kappa = daeParameter("kappa", unit(), self,
                "Gradient penalty")
        self.alpha = daeParameter("alpha", unit(), self,
                " Charge transfer coefficient")
        self.Tabs = daeParameter("Tabs", unit(), self,
                "Temperature in K")
        self.Tref = daeParameter("Tref", unit(), self,
                "Reference temperature in K")
        self.T = daeParameter("T", unit(), self,
                "Non dimensional temperature")
        self.k = daeParameter("k", unit(), self,
                "Boltzmann constant")
        self.currset = daeParameter("C_rate", unit(), self,
                "Discharge C-rate")
        self.cwet = daeParameter("c_wet", unit(), self,
                "Wetted surface concentration")
        self.k0 = daeParameter("k0", unit(), self,
                "exchange current density rate constant")
        self.aO = daeParameter("aO", unit(), self,
                "activity of the oxidized state")
        self.Lx = daeParameter("Lx", unit(), self,
                "Length of particle")
        self.Nx = daeParameter("Nx", unit(), self,
                "Number of volumes in particle")
        self.delx = daeParameter("delx", unit(), self,
                "size of discretization")

    def DeclareEquations(self):
        print "DeclareEquations()"
        daeModel.DeclareEquations(self)

        N = self.N.NumberOfPoints

        # Prepare the noise
        # XXX -- maybe this should be a parameter?
        numnoise = tsteps
        noise_prefac = 1e-3
        noise_data = noise_prefac*np.random.randn(numnoise, N)
        self.noise_vec = np.empty(N, dtype=object)
        self.noise_vec[:] = [noise("Noise", self, unit(), Time(),
                noise_data[:, i], numnoise, self.currset) for i in range(N)]

        # Prepare to evalute the RHS function
        cs = np.empty(N, dtype=object)
        cs[0:N] = [self.c(i) for i in range(N)]
        print "about to set up RHS"
        RHS = self.calc_dcs_dt(cs)

        # The average concentration in the particle
        eq = self.CreateEquation("cbar")
        eq.Residual = self.cbar() - np.sum(cs)/N
        eq.BuildJacobianExpressions = True
        eq.CheckUnitsConsistency = False
        # The following will iterate over rows in the matrix M and
        # create the ∑M(i,:) * y(:) product
        dydt = np.empty(N, dtype=object)
        dydt[0:N] = [self.c.dt(i) for i in range(N)]
#        dydt[N] = self.phi.dt
        Mdydt = self.MX(self.M, dydt)
        # Loop through and make N equations
        print "Forming c Differential Equations"
        for i in range(N):
            # Finally create an equation and set its residual to 'res'
            eq = self.CreateEquation("dydt({num})".format(num=i))
#            # With noise
#            eq.Residual = Mdydt[i] - RHS[i] - self.noise_vec[i]()
            # No noise
            eq.Residual = Mdydt[i] - RHS[i]
            eq.CheckUnitsConsistency = False
        # Total Current Constraint Equation
        # Full equation: self.psivec*dcdt = sum(self.psivec)*self.I
        # parts
        dcdt_vec = self.curr_weight_vec*dydt[0:N]
#        splits_dcdt = self.gensplits(dcdt_vec, num_phi)
#        for i in range(num_phi):
#            dcdt_part = splits_dcdt.next()
#            eq = self.CreateEquation("current_part{num}".format(num=i))
#            eq.Residual = self.phi_parts(i) - self.SUM(dcdt_part)
#            eq.CheckUnitsConsistency = False
        # combined, sum(dc/dt) = I
        eq = self.CreateEquation("Total_Current_Constraint")
        eq.Residual = (
#                Sum(self.phi_parts.array([]), isLargeArray=True)
#                - current_tot)
                np.sum(dcdt_vec) - self.currset())
#                self.SUM(dcdt_vec) - self.currset())
#                - self.I()*np.sum(self.psivec) )
#                Sum(self.c.dt_array([]), isLargeArray=True)
#                - self.I()*np.sum(self.psivec) )
        eq.CheckUnitsConsistency = False

    def calc_dcs_dt(self, cs):
        part_steps = self.N.NumberOfPoints
        cstmp = np.empty(part_steps+2, dtype=object)
        cstmp[1:-1] = cs
        cstmp[0] = self.cwet()
        cstmp[-1] = self.cwet()
        dxs = 1./part_steps
        curv = np.diff(cstmp, 2)/(dxs**2)
#        meanfill = np.sum(cs)/part_steps
        mu = ( self.mu_reg_sln(cs) - self.kappa()*curv
                + self.b()*(cs - self.cbar()) )
        # XXX -- Ask about exp
        act = np.exp(mu)
#        act = self.EXP(mu)
        ecd = ( self.k0() * self.aO()**(1-self.alpha())
                * act**(self.alpha()) * (1-cs) )
        eta = mu - self.phi()
        return ( ecd*(np.exp(-self.alpha()*eta)
                - np.exp((1-self.alpha())*eta)) )
#        return ( ecd*(self.EXP(-self.alpha()*eta)
#                - self.EXP((1-self.alpha())*eta)) )
    def mu_reg_sln(self, c):
        if (type(c[0]) == pyCore.adouble):
            isAdouble = True
        else:
            isAdouble = False
        if isAdouble:
            return np.array([ self.a()*(1-2*c[i])
                    + self.T()*Log(c[i]/(1-c[i]))
                    for i in range(len(c)) ])
        else:
            return ( self.a.GetValue()*(1-2*c)
                    + self.T.GetValue()*np.log(c/(1-c)) )
    def MX(self, mat, objvec):
        if type(mat) is not sprs.csr.csr_matrix:
            raise Exception("MX function designed for csr mult")
        if (type(objvec[0]) == pyCore.adouble):
            isAdouble = True
        else:
            isAdouble = False
        n = objvec.shape[0]
        if isAdouble:
            out = np.empty(n, dtype=object)
        else:
            out = np.zeros(n, dtype=float)
        # Loop through the rows
        for i in range(n):
            low = mat.indptr[i]
            up = mat.indptr[i+1]
            if up > low:
                out[i] = np.sum(
                        mat.data[low:up] * objvec[mat.indices[low:up]] )
            else:
                out[i] = 0.0
        return out
#    def SUM(self, vec):
#        if (type(vec[0]) == pyCore.adouble):
#            return Sum(adouble_array.FromNumpyArray(vec),
#                isLargeArray=True)
#        else:
#            return np.sum(vec)
#    def EXP(self, vec):
#        if (type(vec[0]) == pyCore.adouble):
#            out = np.empty(len(vec), dtype=object)
#            return np.array([Exp(vec[i]) for i in range(len(vec))])
#        else:
#            return np.exp(vec)

class simRay(daeSimulation):
    def __init__(self):
        print "initalize simRay"
        daeSimulation.__init__(self)
        self.m = modRay("Ray")

    def SetUpParametersAndDomains(self):
        print "SetUpParametersAndDomains"
        self.m.Tabs.SetValue(T)
        self.m.Tref.SetValue(Tref)
        self.m.T.SetValue(float(T)/Tref)
        self.m.k.SetValue(k)
        self.m.alpha.SetValue(alpha)
        self.m.a.SetValue(dim_a/(k*Tref))
        self.m.b.SetValue(dim_b/(k*Tref*rhos))
        self.m.kappa.SetValue(dim_kappa/(k*Tref*rhos*part_size**2))
        self.m.currset.SetValue(currset)
        self.m.cwet.SetValue(cwet)
        self.m.k0.SetValue(1.)
        self.m.aO.SetValue(1.)
        part_steps = int(np.ceil(part_size/solid_disc))
        self.m.Nx.SetValue(part_steps)
        self.m.Lx.SetValue(part_size)
        self.m.delx.SetValue(solid_disc)
        self.m.N.CreateArray(part_steps)

        self.m.M = sprs.lil_matrix((part_steps, part_steps))
        self.m.M.setdiag(np.ones(part_steps))
#        self.m.M[part_steps, 0:N] = 1./part_steps
        self.m.curr_weight_vec = (1./part_steps)*np.ones(part_steps)
        self.m.M = self.m.M.tocsr()

    def SetUpVariables(self):
        print "SetUpVariables"
        N = self.m.N.NumberOfPoints
        # Set initial concentration conditions
        for i in range(N):
            self.m.c.SetInitialCondition(i, 0.01)
        # Guess the initial potential
        self.m.phi.SetInitialGuess(0.0)
        print "done SetUpVariables"

class MyMATDataReporter(daeMatlabMATFileDataReporter):
    """
    See Source code for pyDataReporting.daeMatlabMATFileDataReporter
    """
    def WriteDataToFile(self):
        mdict = {}
#        print self.Process
#        print type(self.Process)
#        print dir(self.Process)
        for var in self.Process.Variables:
            mdict[var.Name] = var.Values
            mdict[var.Name + '_times'] = var.TimeValues
#            mdict['psimask'] = self.Process.psimask
        try:
            scipy.io.savemat(self.ConnectionString,
                             mdict,
                             appendmat=False,
                             format='5',
                             long_field_names=False,
                             do_compression=False,
                             oned_as='row')
        except Exception, e:
            print 'Cannot call scipy.io.savemat(); is SciPy installed?\n' + str(e)

def setupDataReporters(simulation):
    """
    Create daeDelegateDataReporter and add data reporters:
     - daeMatlabMATFileDataReporter
    """
    print "setupDataReporter()"
    datareporter = daeDelegateDataReporter()
    simulation.dr = MyMATDataReporter()
    datareporter.AddDataReporter(simulation.dr)
    # Connect data reporters
    simName = simulation.m.Name + strftime(" [%d.%m.%Y %H:%M:%S]", localtime())
#    matfilename = tempfile.gettempdir() + "/daet.mat"
    matfilename = os.path.join(os.getcwd(),
            "acr_sp_{curr}C.mat".format(curr=currset))
    if (simulation.dr.Connect(matfilename, simName) == False):
        sys.exit()
    return datareporter

def consoleRun():
    print "START"
    # Create Log, Solver, DataReporter and Simulation object
    log          = daePythonStdOutLog()
    daesolver    = daeIDAS()
    simulation   = simRay()
    datareporter = setupDataReporters(simulation)

    # Use SuperLU direct sparse LA solver
    print "Set up LU solver"
    lasolver = pySuperLU.daeCreateSuperLUSolver()
    daesolver.SetLASolver(lasolver)
    
    # Enable reporting of all variables
    simulation.m.SetReportingOn(True)

    # Set the time horizon and the reporting interval
    print "set up simulation time parameters"
    simulation.TimeHorizon = 0.98/abs(currset)
    simulation.ReportingInterval = simulation.TimeHorizon/tsteps

    # Connect data reporter
    print "connect to data reporter"
    simName = simulation.m.Name + strftime(" [%d.%m.%Y %H:%M:%S]", localtime())
    if(datareporter.Connect("", simName) == False):
        sys.exit()

    # Initialize the simulation
    print "initialize"
    simulation.Initialize(daesolver, datareporter, log)

    # Solve at time=0 (initialization)
    print "Initialize the system at t = 0"
    simulation.SolveInitial()

    # Run
    print "run the simulation"
    simulation.Run()
    simulation.Finalize()
    
    print 'Number of numpy.interp1d calls in noise ext. functions:'
    print [n.counter for n in simulation.m.noise_vec]

if __name__ == "__main__":
    import timeit
    time_tot = timeit.timeit("consoleRun()",
            setup="from __main__ import consoleRun",  number=1)
    print "Total time:", time_tot, "s"
#    import cProfile
#    cProfile.run("consoleRun()")
#    consoleRun()
