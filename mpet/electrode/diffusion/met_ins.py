import numpy as np
def met_ins(y):
    k = 100
    return (y*(1-y))*(0.1 + 0.9 / (1 + np.exp(-k * (y-0.2))))

