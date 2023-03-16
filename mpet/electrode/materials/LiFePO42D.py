import numpy as np


def LiFePO42D(self, c_mat, ybar, muR_ref):
    muRtheta = -self.eokT*3.422
    Omega = self.get_trode_param("Omega_a")
    kappax = self.get_trode_param("kappa_x")
    kappay = self.get_trode_param("kappa_y")
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

        # curvy = np.empty(Ny, dtype=object)

        curvy = np.empty_like(c_mat_tmp)
        
        curvy[:,0] = 2 * c_mat_tmp[:,1] - 2 * c_mat_tmp[:,0]
        curvy[:,-1] = 2 * c_mat_tmp[:,-2] - 2 * c_mat_tmp[:,-1]
        curvy[:,1:-1] = np.diff(c_mat_tmp, n=2, axis=1)

        curvx = np.diff(c_mat_tmp, n=2, axis=0)
        # muRhomog = self.T * np.log(c_mat_tmp / (1 - c_mat_tmp)) + Omega * (1 - 2 * c_mat)                                                             
        y_vert_avg = np.average(c_mat_tmp, axis=1)
        y_oriz_avg = np.average(c_mat_tmp, axis=0)
        for i in np.arange(1,Nx+1):
            #     curvy[0] = 2*c_mat_tmp[i,1] - 2*c_mat_tmp[i,0]
            #     curvy[-1] = 2*c_mat_tmp[i,-2] - 2*c_mat_tmp[i,-1]
            # # + 2*dys*beta_s*c_mat_tmp[i-1,-1]*(1-c_mat_tmp[i-1,-1])
            #     curvy[1:Ny-1] = np.diff(c_mat_tmp[i,:],2)
            for j in np.arange(Ny):
                y = c_mat_tmp[i,j]
                muRhomog = self.T*np.log(y/(1-y)) + Omega*(1-2*y)
                muR_mat[i-1,j] = muRhomog
                # curvx = (c_mat_tmp[i-1,j] - c_mat_tmp[i,j]) - (c_mat_tmp[i,j] - c_mat_tmp[i+1,j])
                curv_en = -kappax*curvx[i-1,j]/(dxs**2) - kappay*curvy[i,j]/(dys**2)
                muR_mat[i-1,j] += curv_en
                muR_mat[i-1,j] += self.get_trode_param("Bx")*(y - y_oriz_avg[j])
                muR_mat[i-1,j] += self.get_trode_param("By")*(y - y_vert_avg[i])
        # print(muR_mat)
        actR_mat = np.exp(muR_mat/self.T)
        muR_mat += muRtheta + muR_ref
        return muR_mat, actR_mat
