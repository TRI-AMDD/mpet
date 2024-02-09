import numpy as np
from scipy.ndimage import gaussian_filter

def mu(x, d, Lxs,T):
    mu = T*np.log(x / (1-d-x))
    kx = 0
    for Lx in Lxs:
        mu += Lx * ((1-2*x)**kx) * ((1-2*x) - 2 * kx * (x * (1-x)/(1-2*x)))
        kx += 1
    return mu


def LiFePO42D(self, c_mat, ybar, T, muR_ref):
    muRtheta = -self.eokT*3.422
    L_Li = [3.97034583, 0.09673699, 1.11037291, -0.15444768]
    # for allowing muR_ref to be calculated
    if np.size(c_mat) == 1:
        y = c_mat
        # muRhomog = self.reg_sln(y, T, self.get_trode_param("Omega_a"))
        d = 0
        muRhomog = mu(y, d, L_Li,T)
        muRnonHomog = self.general_non_homog(y, ybar)
        muR = muRhomog + muRnonHomog
        actR = np.exp(muR/T)
        muR += muRtheta + muR_ref
        return muR, actR
    else:
        beta_s = self.get_trode_param("beta_s")
        Ny = np.size(c_mat, 1)
        Nx = np.size(c_mat, 0)
        # defect mat matrix Nx x Ny of random numbers between 0 and 0.06
        # mean = 0.05
        # stddevs = 0.03
        # var = stddevs**2
        # mu = np.log(mean**2/np.sqrt(var+mean**2))
        # sigma = np.sqrt(np.log(var/mean**2 + 1))
        # defect_mat = np.random.lognormal(mu, sigma, (Nx,Ny))
        defect_mat = np.random.rand(Nx,Ny)*0
        # be sure it is not smaller than 0
        # gaussian smoothing
        # defect_mat = gaussian_filter(defect_mat, sigma=0.5)
        dys = 1./(Ny-1)
        dxs = 1./Nx
        ywet = self.get_trode_param("cwet")*np.ones((1,Ny), dtype=object)
        muR = np.zeros((Nx,Ny), dtype=object)
        actR_mat = np.zeros((Nx,Ny), dtype=object)
        # vertically the CHR model requires different boundary conditions
        curvy = np.empty_like(c_mat)
        curvy[:,0] = (2 * c_mat[:,1] - 2 * c_mat[:,0])/(dys**2)
        curvy[:,1:Ny-1] = (np.diff(c_mat, n=2, axis=1))/(dys**2)
        curvy[:,Ny-1] = (2 * c_mat[:,-2] - 2 * c_mat[:,-1] + 2*dys*6*beta_s*c_mat[:,-1]*(1-c_mat[:,-1]))/(dys**2)
        # orizontally the ACR is model just require 2 ghost points
        curvx = np.diff(np.concatenate((ywet,c_mat,ywet), axis=0), n=2, axis=0)/(dxs**2)
        
        # regular solution
        # muR = T*np.log(c_mat/(1- defect_mat -c_mat)) + self.get_trode_param("Omega_a")*(1- defect_mat -2*c_mat)
        muR = mu(c_mat, defect_mat, L_Li,T)
        b = 0.3
        pi = 3.14159
        a = -0.11
        phi = -700
        muR_meta_1 = a*20*pi*np.cos(pi*c_mat-b)*(np.sin(pi*c_mat-b))**19
        muR_meta_2 = phi*6*((c_mat*(1-c_mat))**5)*(1-2*c_mat)
        muR += 0*(muR_meta_1 + muR_meta_2)

        # non-homogeneous
        muR += -self.get_trode_param("kappa_x")*curvx - self.get_trode_param("kappa_y")*curvy
        if self.get_trode_param("mechanics") == False:
            y_vert_avg = np.average(c_mat, axis=1)
            y_oriz_avg = np.average(c_mat, axis=0)
            muR += self.get_trode_param("Bx")*np.subtract(c_mat,y_oriz_avg)
            for i in np.arange(Nx):
                y = c_mat[i,:]
                muR[i,:] += self.get_trode_param("By")*(y - y_vert_avg[i])
        actR_mat = np.exp(muR/T)
        muR += muRtheta + muR_ref
        return muR, actR_mat
