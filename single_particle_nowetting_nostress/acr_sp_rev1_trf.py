# Import the necessary modules
import numpy as N
import pylab as P
import pdb as pdb
from assimulo.solvers import IDA
from assimulo.problem import Implicit_Problem

# Some values we need for the computation
steps = 100
dx = 1.0/steps
kappa = 0.001
a = 4.5
cs0 = 0.01
currset = 0.001
io = 0.1

# First a couple functions to make calculations easier
def chempot(cs):
    divcstmp = N.zeros(cs.shape[0]+2)
    divcstmp[1:-1] = cs
    divcstmp[0] = divcstmp[1]   
    divcstmp[-1] = divcstmp[-2] 
    mu = N.log(cs/(1-cs)) + a*(1-2*cs) - kappa*N.diff(divcstmp,2)/(dx**2)
    return mu
    
def residual(t,y,yd):
    # Calculate our reaction rate
    #pdb.set_trace()
    out = N.zeros(y.shape[0])    
    muvec = chempot(y[0:-1])
    act = N.exp(muvec)
    R = -2.0*io*N.sqrt(act)*(1-y[0:-1])*N.sinh((muvec-y[-1])/2.0)
    out[0:-1] = yd[0:-1]-R
    out[-1] = N.sum(R)*dx-currset
    
    return out
    
t0 = 0.0
y0 = cs0*N.ones(steps+1)
y0[-1] = chempot(cs0*N.ones(1))
yd0 = N.zeros(steps+1)

tfinal = .95*(1.0/currset)
ncp = 200

residual(t0,y0,yd0)

model = Implicit_Problem(residual,y0,yd0,t0)
sim = IDA(model)

algvars = N.ones(steps+1)
algvars[-1] = 0.0
sim.algvar = algvars

sim.make_consistent('IDA_YA_YDP_INIT')

try:
	t, y, yd = sim.simulate(tfinal,ncp)
except:
	pdb.set_trace()
	
P.plot(t,-y[:,-1])
P.show()


