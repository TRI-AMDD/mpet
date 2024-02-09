import numpy as np

def mech_LiFePO4_minus():
    # rotatated -45 degrees
    Cij = np.array([
        [152.5, 43, 54.4,   0, -0.85,  0],
        [43, 175.8, 43,     0, -10.32, 0],
        [54.4, 43, 152.5,   0, -0.85,  0],
        [0,     0,   0,     44.7, 0,    6.9],
        [-0.85,-10.32,-0.85,  0,  52.25, 0],
        [0,   0,   0,       6.9,   0,  44.7]
        ])

    strain_tensor = np.array([
        [0.0125,    0,   -0.0375],
        [0,         0.028,     0],
        [-0.0375,   0,    0.0125]
    ])

    e0 = np.array([strain_tensor[0,0], strain_tensor[1,1], strain_tensor[2,2],
                   strain_tensor[1,2], strain_tensor[0,2], strain_tensor[0,1]])

    return Cij*1e9, e0