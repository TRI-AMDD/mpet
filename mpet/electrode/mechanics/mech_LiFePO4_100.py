import numpy as np

def mech_LiFePO4_100():
    # FePo4 elastic constants (GPa)
    c11 = 157.4
    c22 = 175.8
    c33 = 154
    c44 = 37.8
    c55 = 49.05
    c66 = 51.6
    c13 = 51.2
    c12 = 53.35
    c23 = 32.7

    Cij = np.array([
        [c11, c12, c13, 0, 0, 0],
        [c12, c22, c23, 0, 0, 0],
        [c13, c23, c33, 0, 0, 0],
        [0,   0,   0, c44, 0, 0],
        [0,   0,   0, 0, c55, 0],
        [0,   0,   0, 0, 0, c66]
        ])

    # strain
    # e01 = 0.0517
    # e02 = 0.0359
    # e03 = -0.0186
    e01 = 0.05
    e02 = 0.028
    e03 = -0.025

    strain_tensor = np.array([
        [e01, 0, 0],
        [0, e02, 0],
        [0, 0, e03]
    ])

    e0 = np.array([strain_tensor[0,0], strain_tensor[1,1], strain_tensor[2,2],
                   strain_tensor[1,2], strain_tensor[0,2], strain_tensor[0,1]])

    return Cij*1e9, e0