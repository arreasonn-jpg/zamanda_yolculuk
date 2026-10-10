"""Physical Ricci scalar of g_ij = e^{4 phi} gtilde_ij (v2.0-beta1 fix).

R = e^{-4 phi} (Rtilde + R^phi)
where
  R^phi_ij = -2 D_i D_j phi - 2 gtilde_ij D^2 phi
             + 4 d_i phi d_j phi - 4 gtilde_ij |d phi|^2
and all derivatives are conformal-covariant.
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


def ricci_scalar_physical(phi, gt, h):
    """Return physical 3-Ricci scalar grid (N,N,N)."""
    Nx, Ny, Nz = phi.shape
    dphi = _grad(phi, h)
    Hphi = _hess(phi, h)
    Gt, _ = christoffel_from_grid(gt, h)
    dGt = dG_from_grid(Gt, h)
    R = np.zeros_like(phi)
    for ix in range(Nx):
        for iy in range(Ny):
            for iz in range(Nz):
                g = gt[ix, iy, iz]
                ginv = np.linalg.inv(g)
                G = Gt[ix, iy, iz]
                dG = dGt[ix, iy, iz]
                Rij_tilde = ricci_tensor(G, dG)
                Rtilde = float(np.einsum('ij,ij->', ginv, Rij_tilde))

                dp = dphi[ix, iy, iz]
                Hp = Hphi[ix, iy, iz]
                DDphi = Hp - np.einsum('kij,k->ij', G, dp)
                lap = float(np.einsum('ij,ij->', ginv, DDphi))
                dpdp = np.outer(dp, dp)
                dp2 = float(dp @ ginv @ dp)
                Rphi_ij = -2.0 * DDphi - 2.0 * g * lap \
                          + 4.0 * dpdp - 4.0 * g * dp2
                Rphi = float(np.einsum('ij,ij->', ginv, Rphi_ij))
                ph = float(phi[ix, iy, iz])
                R[ix, iy, iz] = np.exp(-4.0 * ph) * (Rtilde + Rphi)
    return R


def hamiltonian_constraint_full(phi, gt, K, At, h, rho_grid=None):
    """H = R + K^2 - A_ij A^ij - 16 pi G rho / c^2 with full R."""
    C2 = (2.99792458e8) ** 2
    G = 6.67430e-11
    R = ricci_scalar_physical(phi, gt, h)
    # A_ij A^ij: use gtilde (Atilde is conformal)
    At_up = np.einsum('...ij,...jk,...kl->...il',
                      np.linalg.inv(gt), At, np.linalg.inv(gt))
    A2 = np.einsum('...ij,...ij->...', At, At_up)
    rho = rho_grid if rho_grid is not None else np.zeros_like(K)
    return R + K**2 - A2 - 16.0 * np.pi * G * rho / C2


__all__ = ["ricci_scalar_physical", "hamiltonian_constraint_full"]
