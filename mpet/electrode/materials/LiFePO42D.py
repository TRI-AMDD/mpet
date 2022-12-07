import numpy as np


def LiFePO42D(self, c_mat, ybar, muR_ref):
    muRtheta = -self.eokT*3.422
    T = 1.
    cbar = ybar
    B = self.get_trode_param("B")
    Omega = self.get_trode_param("Omega_a")
    kappa = self.get_trode_param("kappa")
    # for allowing muR_ref to be calculated
    if np.size(c_mat) == 1:
        c_mat_tmp = c_mat*np.ones((4,4))
        c_mat = c_mat_tmp
    Ny = np.size(c_mat, 0)
    Nx = np.size(c_mat, 1)
    dys = 1./Ny
    dxs = 1./Nx
    ywet = 0.98*np.ones(Ny+2, dtype=object)
    c_mat_tmp = np.zeros((Ny+2,Nx+2), dtype=object)
    muR_mat = np.empty((Ny,Nx), dtype=object)
    actR_mat = np.empty((Ny,Nx), dtype=object)
    c_mat_tmp[2:,1:-1] = c_mat
    c_mat_tmp[:,-1] = ywet
    c_mat_tmp[:,0] = ywet
    c_mat_tmp[0,2:] = c_mat[0,:]
    c_mat_tmp[1,2:] = c_mat[0,:]
    for i in range(2,Ny+2):
        for j in range(1,Nx+1):
            y = c_mat_tmp[i,j]
            muRhomog = T*np.log(y/(1-y)) + Omega*(1-2*y)
            muR_mat[i-2,j-1] = muRhomog
            curvx = (c_mat_tmp[i,j-1] - c_mat_tmp[i,j]) - (c_mat_tmp[i,j] - c_mat_tmp[i,j+1])
            curvy = (c_mat_tmp[i-2,j] - c_mat_tmp[i-1,j]) - (c_mat_tmp[i-1,j] - c_mat_tmp[i,j])
            curv_en = -kappa*curvx/(dxs**2) - kappa*curvy/(dys**2)
            muR_mat[i-2,j-1] += curv_en
            muR_mat[i-2,j-1] += B*(y - cbar)
            actR_mat[i-2,j-1] = np.exp(muR_mat[i-2,j-1]/T)
            muR_mat[i-2,j-1] += muRtheta + muR_ref
    return muR_mat, actR_mat
