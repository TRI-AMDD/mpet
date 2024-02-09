import numpy as np

def tanh(e,c):
    coeff = 3
    eps_a = np.tanh(coeff*(c-0.5))
    eps_0 = np.tanh(coeff*(-0.5))
    eps_1 = np.tanh(coeff*(0.5))
    eps = (eps_a-eps_0)/(eps_1-eps_0)
    return e*eps