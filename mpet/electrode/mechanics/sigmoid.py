import numpy as np

def sigmoid(e, c):
    return e*(1/(1+np.exp(-10*(c-0.5))))
