#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import tempfile

# Import the necessary modules
import numpy as np

from daetools.pyDAE import *
from daetools.pyDAE.data_reporters import *
from time import localtime, strftime
from daetools.solvers.superlu import pySuperLU
from pyUnits import s

steps = 100
currset = 0.001

class modRay(daeModel):
    def __init__(self, Name, Parent=None, Description=""):
        print "Init modRay"
        daeModel.__init__(self, Name, Parent, Description)

        # Domain where y and f are distributed (in this example 3 items)
        self.N = daeDomain("N", self, unit(), "")
        
        # y and f distributed on N
#        # Be careful about the units
        self.cs = daeVariable("y", no_t, self,
                "Solid concentration profile", [self.N])
        self.phi = daeVariable("phi", no_t, self,
                "potential in the solid/conducting phase")

        # variables we'll reference here
        # To do this correctly, use daeParameter as in "Getting
        # started with daetools" in documentation
        self.dx = 1.0/steps
        self.kappa = 0.001
        self.a = 4.5
        self.cs0 = 0.01
#        self.currset = 0.001
        self.io = 0.1
        self.pflag = False

    def DeclareEquations(self):
        print "DeclareEquations()"
        daeModel.DeclareEquations(self)

        n = self.N.NumberOfPoints

        # Your function f(y,t) s.t. M*(dy/dt) = f(y,t)
        cs = np.empty(n, dtype=object)
        cs[:] = [self.cs(i) for i in range(n)]
        print "setting up RHS:"
        muvec = self.chempot(cs)
        muvec_adb = adouble_array.FromNumpyArray(muvec)
        act_adb = Exp(muvec_adb)
#        act = np.empty(n, dtype=object)
#        act[:] = [act_adb(i) for i in range(n)]
        #
#        cs_adouble = Array([self.y(i) for i in range(n)])
#        cs_adb = adouble_array.FromNumpyArray(cs)
#        phi = self.(n)
#        print muvec_adb-self.phi()
#        R = ( -2.0*self.io*Sqrt(act_adb)*(1-self.cs.array([]))
#                    *Sinh((muvec_adb-self.phi())/2.0) )
        R = ( -self.io*Sqrt(act_adb)*(1-self.cs.array([]))
                    *(Exp((muvec_adb-self.phi())/2.0)
                    - Exp(-(muvec_adb-self.phi())/2.0))
                    )
        # Loop through and make N equations
        print "Forming Equations"
        for i in range(n):
#            res = Constant(0 / s) # For unit consistency only
#            res = Mdydt[i] - RHSofy[i]# * Constant(1 / s)
#            R = ( -2.0*self.io*Sqrt(act(i))*(1-y(i))
#                    *Sinh((muvec(i)-y(n))/2.0) )
            # Finally create an equation and set its residual to 'res'
            eq = self.CreateEquation("dydt({num})".format(num=i))
#            print "made dydt({num}) of {tot}".format(num=i,tot=n)
            eq.Residual = self.cs.dt(i) - R(i)
            eq.CheckUnitsConsistency = False
        eq = self.CreateEquation("current_integral_constraint")
        eq.Residual = Sum(R)*self.dx - currset
        eq.CheckUnitsConsistency = False
        print "Done with Equations"

    def chempot(self, cs):
        n = cs.shape[0]
        divcstmp = np.empty(n + 2, dtype=object)
#        divcstmp[1:-1] = [cs[i] for i in range(n)]
        divcstmp[1:-1] = cs
        divcstmp[0] = divcstmp[1]   
        divcstmp[-1] = divcstmp[-2] 
#        mu = (np.log(cs/(1-cs)) + self.a*(1-2*cs)
#                - self.kappa*np.diff(divcstmp,2)/(self.dx**2))
        grad_term = self.kappa*np.diff(divcstmp, 2)/(self.dx**2)
        mu = np.array([
            Log(cs[i]/(1-cs[i])) + self.a*(1-2*cs[i]) - grad_term[i] 
            for i in range(n) ])
        return mu

class simRay(daeSimulation):
    def __init__(self):
        print "initalize simRay"
        daeSimulation.__init__(self)
        self.m = modRay("Ray")
#        self.steps = 100
#        self.n = self.steps+1

    def SetUpParametersAndDomains(self):
        print "SetUpParametersAndDomains"
        self.m.N.CreateArray(steps)

    def SetUpVariables(self):
        print "SetUpVariables"
        # Initial conditions
        for i in range(steps):
            self.m.cs.SetInitialCondition(i, 0.01)
        self.m.phi.SetInitialGuess(0.0)

class MyDataReporter(daeDataReporterLocal):
    """
    See Tutorial 8
    """
    def __init__(self):
        print "init MyDataReporter"
        daeDataReporterLocal.__init__(self)
        self.ProcessName=""
    def Connect(self, ConnectionString, ProcessName):
        print "Connect()"
        self.ProcessName = ProcessName
        try:
            self.f = open(ConnectionString, "w")
        except IOError:
            return False
        return True
    def Disconnect(self):
        print "Disconnect()"
        self.Write()
        return True
    def MakeString(self):
        print "MakeString()"
        s = "Process name: " + self.ProcessName + "\n"
        variables = self.Process.Variables
        for var in variables:
            values = var.Values
            domains = var.Domains
            times = var.TimeValues
            s += " - Variable: " + var.Name + "\n"
            s += "    - Domains:" + "\n"
            for domain in domains:
                s += "       - " + domain.Name + "\n"
            s += "    - Values:" + "\n"
            for i in range(len(times)):
                s += "      - Time: " + str(times[i]) + "\n"
                s += "        " + str(values[i, ...]) + "\n"
        return s
    def Write(self):
        print "Write()"
        try:
            content = self.MakeString()
#            print content
            self.f.write(content)
            self.f.close()
        except IOError:
            self.f.close()
            return False
    def IsConnected(self):
        print "IsConnected()"
        return True

def setupDataReporters(simulation):
    """
    Create daeDelegateDataReporter and add data reporters:
     - MyDataReporterLocal
     - daeTCPIPDataReporter
     - daeMatlabMATFileDataReporter
     - daePlotDataReporter
    """
    print "setupDataReporter()"
    datareporter = daeDelegateDataReporter()
    simulation.dr1 = MyDataReporter()
#    simulation.dr2 = daeTCPIPDataReporter()
    simulation.dr3 = daeMatlabMATFileDataReporter()
    simulation.dr4 = daePlotDataReporter()
    datareporter.AddDataReporter(simulation.dr1)
#    datareporter.AddDataReporter(simulation.dr2)
    datareporter.AddDataReporter(simulation.dr3)
    datareporter.AddDataReporter(simulation.dr4)
    # Connect data reporters
    simName = simulation.m.Name + strftime(" [%d.%m.%Y %H:%M:%S]", localtime())
    filename = tempfile.gettempdir() + "/acr_sp.out"
    matfilename = tempfile.gettempdir() + "/acr_sp.mat"
    txtfilename = tempfile.gettempdir() + "/acr_sp.txt"
    if (simulation.dr1.Connect(filename, simName) == False):
        sys.exit()
#    if (simulation.dr2.Connect("", simName) == False):
#        sys.exit()
    if (simulation.dr3.Connect(matfilename, simName) == False):
        sys.exit()
    if (simulation.dr4.Connect(txtfilename, simName) == False):
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
    simulation.TimeHorizon = 0.95*(1.0/currset)
    simulation.ReportingInterval = simulation.TimeHorizon/100

    # Connect data reporter
    print "connect to data reporter"
    simName = simulation.m.Name + strftime(" [%d.%m.%Y %H:%M:%S]", localtime())
    if(datareporter.Connect("", simName) == False):
        sys.exit()

    # Initialize the simulation
    print "initialize"
    simulation.Initialize(daesolver, datareporter, log)

    # Save the model report and the runtime model report
    print "set to save report and runtime model"
#    simulation.m.SaveModelReport(simulation.m.Name + ".xml")
    #print "Save runtime model report"
#    simulation.m.SaveRuntimeModelReport(simulation.m.Name + "-rt.xml")

    # Solve at time=0 (initialization)
    print "Initialize the system at t = 0"
    simulation.SolveInitial()

    # Run
    print "run the simulation"
    simulation.Run()
    simulation.Finalize()

if __name__ == "__main__":
    consoleRun()
