#!/usr/bin/env python

"""'Hello, world!' example."""

import os
import sys
from time import localtime, strftime

import daetools.pyDAE as dae
from daetools.pyDAE.data_reporters import daeMatlabMATFileDataReporter
import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio


class modTutorial(dae.daeModel):
    """Define modTutorial class."""

    def __init__(self, Name, Parent=None, Description=''):
        """Call __init__ of base class and define variables."""
        dae.daeModel.__init__(self, Name, Parent, Description)
        self.tau = dae.daeVariable('tau', dae.no_t, self, 'Time elapsed in the process')

    def DeclareEquations(self):
        """Call DeclareEquations of base class and define equations."""
        dae.daeModel.DeclareEquations(self)
        eq = self.CreateEquation('Time', 'Differential equation to calculate the time elapsed in the process')
        eq.Residual = self.tau.dt() - 1.0
        eq.CheckUnitsConsistency = False


class simTutorial(dae.daeSimulation):
    """Define simTutorial class."""

    def __init__(self):
        """Call __init__ of base class and define models."""
        dae.daeSimulation.__init__(self)
        self.m = modTutorial('hello_world')

    def SetUpParametersAndDomains(self):
        """Set up parameters and domains."""
        pass

    def SetUpVariables(self):
        """Set up variables."""
        self.m.tau.SetInitialCondition(0.0)


def setupDataReporters(simulation):
    """Define setupDataReporters function."""
    datareporter = dae.daeDelegateDataReporter()
    simulation.dr = daeMatlabMATFileDataReporter()
    datareporter.AddDataReporter(simulation.dr)
    simName = simulation.m.Name + strftime(' [%d.%m.%Y %H:%M:%S]', localtime())
    matfilename = os.getcwd() + '/' + simulation.m.Name + '.mat'
    if not simulation.dr.Connect(matfilename, simName):
        sys.exit()
    return datareporter


def consoleRun(config):
    """Initialize and run simulations."""
    log = dae.daePythonStdOutLog()
    daesolver = dae.daeIDAS()
    simulation = simTutorial()
    datareporter = setupDataReporters(simulation)
    simulation.TimeHorizon = config['t_final']
    simulation.ReportingInterval = config['t_interval']
    simulation.m.SetReportingOn(True)

    simulation.Initialize(daesolver, datareporter, log)
    simulation.SolveInitial()
    simulation.Run()
    simulation.Finalize()

if __name__ == '__main__':
    N_t = 11
    t_final = 10

    config = {}
    config['t_final'] = t_final
    config['t_interval'] = t_final / (N_t-1)
    consoleRun(config)

    # Plot results
    data = sio.loadmat('hello_world')
    os.remove('hello_world.mat')
    tau = data['hello_world.tau']
    t = np.linspace(0, t_final, N_t)
    tau = tau.squeeze()
    plt.plot(t, tau)
    plt.xlabel(r'$t$')
    plt.ylabel(r'$\tau$')
    plt.axis('tight')
    plt.show()
