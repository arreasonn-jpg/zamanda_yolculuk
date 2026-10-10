"""Full BSSN RHS with covariant Hessian + D^2 alpha (v2.0-beta1).

Reference: Alcubierre 2008, "Introduction to 3+1 NR", §11.3.
Shift beta = 0 (to be added in beta2).

Key corrections vs v1.x:
- dK now contains -D^2 alpha (gauge back-reaction).
- dAtilde uses covariant Hessian: D_i D_j phi = H_ij - Gamma^k_ij d_k phi.
- dAtilde trace-free part subtracted w.r.t. conformal metric.
"""
from __future__ import annotations
import numpy as np
from .bssn_grid import d_central, christoffel_from_grid, dG_from_grid
from .bssn_rhs import ricci_tensor


def _grad(f, h):
    return np.stack(d_central(f, h), axis=-1)


def _hess(f, h):
    d = d_central(f, h)
    H = np.zeros(f.shape + (3, 3))
    for i in range(3):
        dd = d_central(d[i], h)
        for j in range(3):
            H[..., i, j] = dd[j]
    return H


def bssn_rhs_beta1(phi, gt, K, At, alpha, h):
    """Full BSSN RHS (shift=0). All fields are (N,N,N) grids."""
    Nx, Ny, Nz = phi.shape
    dphi = _grad(phi, h)
    Hphi = _hess(phi, h)
    Halpha = _hess(alpha, h)
    dalpha = _grad(alpha, h)
    dK_grad = _grad(K, h)
    G_tilde, _ = christoffel_from_grid(gt, h)
    dG_tilde = dG_from_grid(G_tilde, h)

    dphi_t = np.zeros_like(phi)
    dgt_t = np.zeros_like(gt)
    dK_t = np.zeros_like(K)
    dAt_t = np.zeros_like(At)
    dG_t = np.zeros(phi.shape + (3,))

    for ix in range(Nx):
        for iy in range(Ny):
            for iz in range(Nz):
                a = alpha[ix,iy,iz]
                Kp = K[ix,iy,iz]
                ph = phi[ix,iy,iz]
                g = gt[ix,iy,iz]
                A = At[ix,iy,iz]
                G = G_tilde[ix,iy,iz]
                dG = dG_tilde[ix,iy,iz]
                ginv = np.linalg.inv(g)
                em4phi = np.exp(-4*ph)

                dphi_p = dphi[ix,iy,iz]
                Hphi_p = Hphi[ix,iy,iz]
                Ha_p = Halpha[ix,iy,iz]
                da_p = dalpha[ix,iy,iz]
                dK_p = dK_grad[ix,iy,iz]

                A_up = ginv @ A @ ginv          # A^ij
                A2_scalar = float(np.sum(A * A_up))
                A2_tensor = A @ ginv @ A        # A_ik A^k_j

                Rij = ricci_tensor(G, dG)

                dphi_t[ix,iy,iz] = -a * Kp / 6.0
                dgt_t[ix,iy,iz]  = -2.0 * a * A

                # dK = -D^2 alpha + alpha (A_ij A^ij + K^2/3)
                DD_alpha = (np.einsum('ij,ij->', ginv, Ha_p)
                            - np.einsum('ij,kij,k->', ginv, G, da_p))
                dK_t[ix,iy,iz] = -DD_alpha + a * (A2_scalar + Kp**2 / 3.0)

                # dAtilde
                DD_phi = Hphi_p - np.einsum('kij,k->ij', G, dphi_p)
                bracket = -2.0*a*DD_phi + 4.0*a*np.outer(dphi_p, dphi_p) + a*Rij
                bracket -= (np.trace(ginv @ bracket) / 3.0) * g
                dAt_val = em4phi * bracket + a * (Kp * A - 2.0 * A2_tensor)
                dAt_val -= (np.trace(ginv @ dAt_val) / 3.0) * g
                dAt_t[ix,iy,iz] = dAt_val

                # dGamma
                t1 = -2.0 * (A_up @ da_p)
                ta = np.einsum('ijk,jk->i', G, A_up)
                tb = -(2.0/3.0) * (ginv @ dK_p)
                tc = 6.0 * (A_up @ dphi_p)
                dG_t[ix,iy,iz] = t1 + 2.0*a*(ta + tb + tc)

    return {"dphi": dphi_t, "dgt": dgt_t, "dK": dK_t,
            "dAt": dAt_t, "dGamma": dG_t}


__all__ = ["bssn_rhs_beta1"]
