"""Helper functions to get information about the mesh/geometry of the simulated particles."""
import numpy as np

import mpet.utils as utils


def get_unit_solid_discr(Shape, N):
    if N == 1:  # homog particle, hopefully
        r_vec = None
        volfrac_vec = np.ones(1)
    elif Shape == "C3":
        r_vec = None
        # For 1D particle, the vol fracs are simply related to the
        # length discretization
        volfrac_vec = (1./N) * np.ones(N)  # scaled to 1D particle volume
    elif Shape == "sphere":
        Rs = 1.  # (non-dimensionalized by itself)
        dr = Rs/(N - 1)
        r_vec = np.linspace(0, Rs, N)
        vol_vec = 4*np.pi*(r_vec**2 * dr + (1./12)*dr**3)
        vol_vec[0] = 4*np.pi*(1./24)*dr**3
        vol_vec[-1] = (4./3)*np.pi*(Rs**3 - (Rs - dr/2.)**3)
        Vp = 4./3.*np.pi*Rs**3
        volfrac_vec = vol_vec/Vp
    elif Shape == "cylinder":
        Rs = 1.  # (non-dimensionalized by itself)
        h = 1.
        dr = Rs / (N - 1)
        r_vec = np.linspace(0, Rs, N)
        vol_vec = np.pi * h * 2 * r_vec * dr
        vol_vec[0] = np.pi * h * dr**2 / 4.
        vol_vec[-1] = np.pi * h * (Rs * dr - dr**2 / 4.)
        Vp = np.pi * Rs**2 * h
        volfrac_vec = vol_vec / Vp
    else:
        raise NotImplementedError("Fix shape volumes!")
    return r_vec, volfrac_vec


def get_dr_edges(shape, N):
    r_vec = get_unit_solid_discr(shape, N)[0]
    dr = edges = None
    if r_vec is not None:
        Rs = 1.
        dr = r_vec[1] - r_vec[0]
        edges = np.hstack((0, utils.mean_linear(r_vec), Rs))
    return dr, edges


def calc_curv(c, dr, r_vec, Rs, beta_s, particleShape):
    N = len(c)
    curv = np.empty(N, dtype=c.dtype)
    # Here, beta_s = n*grad(c) at the surface.
    # beta_s = (c_N - c_{N-2})/(2*dr)
    # c_N = c_{N-2} + 2*dr*beta_s
    if particleShape == "sphere":
        curv[0] = 3 * (2*c[1] - 2*c[0]) / dr**2
        curv[1:N-1] = (
            np.diff(c, 2)/dr**2
            + (2./r_vec[1:-1])*(c[2:] - c[0:-2])/(2*dr))
        curv[N-1] = (
            (2./Rs)*beta_s
            + (2*c[-2] - 2*c[-1] + 2*dr*beta_s)/dr**2)
    elif particleShape == "cylinder":
        curv[0] = 2 * (2*c[1] - 2*c[0]) / dr**2
        curv[1:N-1] = (
            np.diff(c, 2)/dr**2
            + (1./r_vec[1:-1])*(c[2:] - c[0:-2])/(2*dr))
        curv[N-1] = (
            (1./Rs)*beta_s
            + (2*c[-2] - 2*c[-1] + 2*dr*beta_s)/dr**2)
#    elif particleShape == "rod":
#        curv[0] = (2*c[1] - 2*c[0]) / dr**2
#        curv[1:N-1] = np.diff(c, 2)/dr**2
#        curv[N-1] = (2*c[-2] - 2*c[-1] + 2*dr*beta_s)/dr**2
    else:
        raise NotImplementedError("calc_curv_c only for sphere and cylinder")
    return curv


def get_elyte_disc(Nvol, L, poros, BruggExp):
    out = {}
    # Width of each cell
    out["dxvec"] = utils.get_dxvec(L, Nvol)

    # Distance between cell centers
    dxtmp = np.hstack((out["dxvec"][0], out["dxvec"], out["dxvec"][-1]))
    out["dxd1"] = utils.mean_linear(dxtmp)

    # The porosity vector
    out["porosvec"] = utils.get_asc_vec(poros, Nvol)
    porosvec_pad = utils.pad_vec(out["porosvec"])

    # Vector of Bruggeman exponents
    Brugg_pad = utils.pad_vec(utils.get_asc_vec(BruggExp, Nvol))

    # Vector of posority/tortuosity (assuming Bruggeman)
    out["eps_o_tau"] = porosvec_pad/porosvec_pad**(Brugg_pad)

    return out
