import numpy as np


def LiFePO42D(self, c_mat, ybar, muR_ref):
    muRtheta = -self.eokT*3.422
    T = 1.
    cbar = ybar
    B = self.get_trode_param("B")
    Omega = self.get_trode_param("Omega_a")
    kappax = self.get_trode_param("kappa")
    kappay = kappax*10
    # beta_s = self.get_trode_param("beta_s")
    # for allowing muR_ref to be calculated
    if np.size(c_mat) == 1:
        y = c_mat
        muRhomog = self.reg_sln(y, self.get_trode_param("Omega_a"))
        muRnonHomog = self.general_non_homog(y, ybar)
        muR = muRhomog + muRnonHomog
        actR = np.exp(muR/self.T)
        muR += muRtheta + muR_ref
        return muR, actR
    else:
        Ny = np.size(c_mat, 1)
        Nx = np.size(c_mat, 0)
        dys = 1./Ny
        dxs = 1./Nx
        ywet = 0.98*np.ones(Ny, dtype=object)
        c_mat_tmp = np.zeros((Nx+2,Ny), dtype=object)
        muR_mat = np.zeros((Nx,Ny), dtype=object)
        actR_mat = np.zeros((Nx,Ny), dtype=object)
        # from second to second to last row we have the original conc
        c_mat_tmp[1:-1,:] = c_mat
        # first and last row is 0.98
        c_mat_tmp[-1,:] = ywet
        c_mat_tmp[0,:] = ywet
        curvy = np.empty(Ny, dtype=object)
        # print(muR_mat)
        # plan: produce the matrix muR without curvy
        # then another for loop for the curvy
        # if then works try to optmize the loops
        for i in range(1,Nx+1):
            curvy[0] = 2*c_mat_tmp[i,1] - 2*c_mat_tmp[i,0]
            curvy[-1] = 2*c_mat_tmp[i,-2] - 2*c_mat_tmp[i,-1] \
                # + 2*dys*beta_s*c_mat_tmp[i-1,-1]*(1-c_mat_tmp[i-1,-1])
            curvy[1:Ny-1] = np.diff(c_mat_tmp[i,:],2)
            for j in range(Ny):
                y = c_mat_tmp[i,j]
                muRhomog = T*np.log(y/(1-y)) + Omega*(1-2*y)
                muR_mat[i-1,j] = muRhomog
                curvx = (c_mat_tmp[i-1,j] - c_mat_tmp[i,j]) - (c_mat_tmp[i,j] - c_mat_tmp[i+1,j])
                # curvy[i-1,j] = (c_mat_tmp[i-1,j-1] - c_mat_tmp[i-1,j]) - (c_mat_tmp[i-1,j]
                #                                                           - c_mat_tmp[i-1,j+1])
                # curvy[i-1,-1] = 2*c_mat_tmp[i-1,-2] - 2*c_mat_tmp[i-1,-1]
                # + 2*dys*beta_s*c_mat_tmp[i-1,-1]*(1-c_mat_tmp[i-1,-1])
                # curv_en = -kappax*curvx/(dxs**2) - kappay*curvy/(dys**2)
                curv_en = -kappax*curvx/(dxs**2) - kappay*curvy[j]/(dys**2)
                muR_mat[i-1,j] += curv_en
                muR_mat[i-1,j] += B*(y - cbar)
        # print(muR_mat)
        actR_mat = np.exp(muR_mat/T)
        muR_mat += muRtheta + muR_ref
        return muR_mat, actR_mat
