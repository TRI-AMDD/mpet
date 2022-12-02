import numpy as np


def LiFePO42D(self, c_mat, ybar, muR_ref):
    muRtheta = -self.eokT*3.422

    Ny = np.size(c_mat, 0)
    dys = 1./Ny
    Nx = np.size(c_mat, 1)
    dxs = 1./Nx
    muR_mat = np.empty((Ny,Nx), dtype=object)
    actR_mat = np.empty((Ny,Nx), dtype=object)
    # shape = self.get_trode_param("shape")
    kappa = self.get_trode_param("kappa")
    B = self.get_trode_param("B")
    # cwet = self.get_trode_param("cwet")
    # ytmp = np.empty(N+2, dtype=object)
    # muRnonHomog_mat = self.general_non_homog(c_mat, ybar)
    for i in range(Ny-2):
        for j in range(Nx-2):
            y = c_mat[i,j]
            muRhomog = self.reg_sln(y, self.get_trode_param("Omega_a"))
            muR_mat[i,j] = muRhomog
            curvx = (c_mat[i,j] - c_mat[i,j+1]) - (c_mat[i,j+1] - c_mat[i,j+2])
            curvy = (c_mat[i,j] - c_mat[i+1,j]) - (c_mat[i+1,j] - c_mat[i+2,j])
            curv_en = -kappa*curvx/(dxs**2) - kappa*curvy/(dys**2)
            muR_mat[i,j] += curv_en
            muR_mat[i,j] += B*(y - ybar)
            actR_mat[i,j] = np.exp(muR_mat[i,j]/self.T)
            muR_mat[i,j] += muRtheta + muR_ref
    return muR_mat, actR_mat
