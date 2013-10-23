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
partsize = 1000                       # particle size, nm

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
        lowerBound=0, upperBound=1, initialGuess=0.25,
        absTolerance=1e-6)
elec_pot_t = daeVariableType(name="elec_pot_t", units=unit(),
        lowerBound=-1e20, upperBound=1e20, initialGuess=0,
        absTolerance=1e-5)

class noise(daeScalarExternalFunction):
    def __init__(self, Name, Model, units, time, time_vec,
            noise_data, previous_output, position):
        arguments = {}
        self.counter = 0
        self.saved = 0
        self.previous_output = previous_output
        self.time_vec = time_vec
        self.noise_data = noise_data
        self.tlo = time_vec[0]
        self.thi = time_vec[-1]
        self.numnoise = len(time_vec)
        arguments["time"] = time
        self.position = position
        daeScalarExternalFunction.__init__(self, Name, Model, units, arguments)

    def Calculate(self, values):
        time = values["time"]
        # A derivative for Jacobian is requested - return always 0.0
        if time.Derivative != 0:
            return adouble(0)
        # Store the previous time value to prevent excessive
        # interpolation.
        if len(self.previous_output) > 0 and self.previous_output[0] == time.Value:
            self.saved += 1
            return adouble(float(self.previous_output[1][self.position]))
        indx = (float(time.Value - self.tlo)/(self.thi-self.tlo) *
                (self.numnoise - 1))
        ilo = np.floor(indx)
        ihi = np.ceil(indx)
        # If we're exactly at a time in time_vec
        if ilo == ihi:
            noise_vec = self.noise_data[ilo, :]
        else:
            noise_vec = (self.noise_data[ilo, :] +
                    (time.Value - self.time_vec[ilo]) /
                    (self.time_vec[ihi] - self.time_vec[ilo]) *
                    (self.noise_data[ihi, :] - self.noise_data[ilo, :])
                    )
#        print 'At time = %f the function with index = %d interpolated the noise' % (time.Value, int(self.position)) 
        # previous_output is a reference to a common object and must
        # be updated here - not deleted.  using self.previous_output = []
        # it will delete the common object and create a new one
        self.previous_output[:] = [time.Value, noise_vec] # it is a list now not a tuple
        self.counter += 1
        return adouble(float(noise_vec[self.position]))

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
        numnoise = tsteps/10
        noise_prefac = 1e-3
        noise_data = noise_prefac*np.random.randn(numnoise, N)
        # a vector going from 0 to the max simulation time.
        time_vec = np.linspace(0, (1./self.currset.GetValue()), numnoise)
        # Previous_output is common for all external functions
        previous_output = []
        # daeScalarExternalFunction (noise interpolation done as vector)
        self.noise_local = np.empty(N, dtype=object)
        self.noise_local[:] = [noise("Noise", self, unit(), Time(),
                                     time_vec, noise_data, previous_output, _position_)
                               for _position_ in range(N)]

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

        # Make a numpy array of the derivative variables for easy
        # manipulation.
        dydt = np.empty(N, dtype=object)
        dydt[0:N] = [self.c.dt(i) for i in range(N)]
        # Loop through and make N equations
        print "Forming c Differential Equations"
        for i in range(N):
            # Finally create an equation and set its residual to 'res'
            eq = self.CreateEquation("dydt({num})".format(num=i))
#            # With noise
#            # scalar function, referencing common noise-vector
#            eq.Residual = (dydt[i] - RHS[i] - self.noise_local[i]())
            # No noise
            eq.Residual = dydt[i] - RHS[i]
            eq.CheckUnitsConsistency = False
        # Total Current Constraint Equation
        # Full equation: self.psivec*dcdt = sum(self.psivec)*self.I
        # parts
        curr_cons_weight_vec = (1./N)*np.ones(N)
        dcdt_vec = curr_cons_weight_vec*dydt[0:N]
        # combined, sum(dc/dt) = I
        eq = self.CreateEquation("Total_Current_Constraint")
        eq.Residual = (
                np.sum(dcdt_vec) - self.currset())
        eq.CheckUnitsConsistency = False

    def calc_dcs_dt(self, cs):
        part_steps = self.N.NumberOfPoints
        cstmp = np.empty(part_steps+2, dtype=object)
        cstmp[1:-1] = cs
        cstmp[0] = self.cwet()
        cstmp[-1] = self.cwet()
        dxs = 1./part_steps
        curv = np.diff(cstmp, 2)/(dxs**2)
        mu = ( self.mu_reg_sln(cs) - self.kappa()*curv
                + self.b()*(cs - self.cbar()) )
        act = np.exp(mu)
        ecd = ( self.k0() * self.aO()**(1-self.alpha())
                * act**(self.alpha()) * (1-cs) )
        eta = mu - self.phi()
        return ( ecd*(np.exp(-self.alpha()*eta)
                - np.exp((1-self.alpha())*eta)) )
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
        for var in self.Process.Variables:
            mdict[var.Name] = var.Values
            mdict[var.Name + '_times'] = var.TimeValues
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
#    lasolver = pySuperLU.daeCreateSuperLUSolver()
    lasolver = pyTrilinos.daeCreateTrilinosSolver("Amesos_Umfpack", "")
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
    
#    print 'Number of numpy.interp1d calls in noise ext. functions:'
#    print simulation.m.noise_local.counter
#    print simulation.m.noise_local.saved
#    print [n.counter for n in simulation.m.noise_vec]

if __name__ == "__main__":
    import timeit
    time_tot = timeit.timeit("consoleRun()",
            setup="from __main__ import consoleRun",  number=1)
    print "Total time:", time_tot, "s"
#    import cProfile
#    cProfile.run("consoleRun()")
#    consoleRun()
