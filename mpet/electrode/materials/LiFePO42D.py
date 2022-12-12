import numpy as np


def LiFePO42D(self, c_mat, ybar, muR_ref):
    muRtheta = -self.eokT*3.422
    Omega = self.get_trode_param("Omega_a")
    kappax = self.get_trode_param("kappa")
    kappay = kappax*10
    # beta_s = self.get_trode_param("beta_s")
    # for allowing muR_ref to be calculated
    # FePo4 elastic constants (GPa)
    c11 = 175.9
    c22 = 153.6
    c33 = 135.0
    c44 = 38.8
    c55 = 47.5
    c66 = 55.6
    c13 = 54.0
    c12 = 29.6
    c23 = 19.6

    Cij = np.zeros((6,6))
    Cij[0,0] = c11
    Cij[1,1] = c22
    Cij[2,2] = c33
    Cij[3,3] = c44
    Cij[4,4] = c55
    Cij[5,5] = c66
    Cij[1,0] = c12
    Cij[0,1] = c12
    Cij[2,0] = c13
    Cij[0,2] = c13
    Cij[1,2] = c23
    Cij[2,1] = c23
    # strain
    e1 = 0.0517
    e2 = 0.0359
    e3 = -0.0186
    ej = np.array([e1, e2, e3, 0, 0, 0])
    sigma_i = np.dot(Cij,ej)
    print(sigma_i)

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
            curvy[-1] = 2*c_mat_tmp[i,-2] - 2*c_mat_tmp[i,-1]
        # + 2*dys*beta_s*c_mat_tmp[i-1,-1]*(1-c_mat_tmp[i-1,-1])
            curvy[1:Ny-1] = np.diff(c_mat_tmp[i,:],2)
            for j in range(Ny):
                y = c_mat_tmp[i,j]
                muRhomog = self.T*np.log(y/(1-y)) + Omega*(1-2*y)
                muR_mat[i-1,j] = muRhomog
                curvx = (c_mat_tmp[i-1,j] - c_mat_tmp[i,j]) - (c_mat_tmp[i,j] - c_mat_tmp[i+1,j])
                curv_en = -kappax*curvx/(dxs**2) - kappay*curvy[j]/(dys**2)
                muR_mat[i-1,j] += curv_en
                muR_mat[i-1,j] += self.get_trode_param("B")*(y - ybar)
        # print(muR_mat)
        actR_mat = np.exp(muR_mat/self.T)
        muR_mat += muRtheta + muR_ref
        return muR_mat, actR_mat
