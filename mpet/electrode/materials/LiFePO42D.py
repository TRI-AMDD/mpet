import numpy as np


def LiFePO42D(self, c_mat, ybar, muR_ref):
    muRtheta = -self.eokT*3.422
    Omega = self.get_trode_param("Omega_a")
    kappax = self.get_trode_param("kappa_x")
    kappay = self.get_trode_param("kappa_y")
    ywet = self.get_trode_param("cwet")
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
        ywet = ywet*np.ones((1,Ny), dtype=object)
        muR_mat = np.zeros((Nx,Ny), dtype=object)
        actR_mat = np.zeros((Nx,Ny), dtype=object)
        curvy = np.empty_like(c_mat)
        curvy[:,0] = (2 * c_mat[:,1] - 2 * c_mat[:,0])/(dys**2)
        curvy[:,1:Ny-1] = (np.diff(c_mat, n=2, axis=1))/(dys**2)
        curvy[:,Ny-1] = (2 * c_mat[:,-2] - 2 * c_mat[:,-1])/(dys**2)
        curvx = np.diff(np.concatenate((ywet,c_mat,ywet), axis=0), n=2, axis=0)/(dxs**2)
        y_vert_avg = np.average(c_mat, axis=1)
        y_oriz_avg = np.average(c_mat, axis=0)
        muR_mat = self.T*np.log(c_mat/(1-c_mat)) + Omega*(1-2*c_mat)
        muR_mat += -kappax*curvx - kappay*curvy
        muR_mat += self.get_trode_param("Bx")*np.subtract(c_mat,y_oriz_avg)
        for i in np.arange(Nx):
            y = c_mat[i,:]
            muR_mat[i,:] += self.get_trode_param("By")*(y - y_vert_avg[i])
        actR_mat = np.exp(muR_mat/self.T)
        muR_mat += muRtheta + muR_ref
        return muR_mat, actR_mat
