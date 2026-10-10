"""Full BSSN evolution of the conformal connection Gamma^i (v2.0-alpha2)."""
from __future__ import annotations
import numpy as np
from .bssn_grid import d_central, christoffel_from_grid


def _grad(f, h):
    """f: (N,N,N) -> (N,N,N,3) array of partials."""
    return np.stack(d_central(f, h), axis=-1)


def _grad3(A, h):
    """A: (N,N,N,3,3) -> (N,N,N,3,3,3), dA[...,i,j,k] = d_k A_ij."""
    Nx, Ny, Nz = A.shape[:3]
    dA = np.zeros((Nx,Ny,Nz,3,3,3))
    for i in range(3):
        for j in range(3):
            comp = A[..., i, j]
            dx, dy, dz = d_central(comp, h)
            dA[..., i, j, 0] = dx
            dA[..., i, j, 1] = dy
            dA[..., i, j, 2] = dz
    return dA


def dGamma_RHS(Gamma, gt, At, K, phi, alpha, h):
    """d_t Gamma^i = -2 Atilde^ij d_j alpha
                    + 2 alpha [ G^i_jk Atilde^jk - (2/3) gt^ij d_j K
                                + 6 Atilde^ij d_j phi ]."""
    dalpha = _grad(alpha, h)     # (N,N,N,3)
    dK = _grad(K, h)
    dphi = _grad(phi, h)
    G, _ = christoffel_from_grid(gt, h)   # (N,N,N,3,3,3)

    # -2 Atilde^ij d_j alpha
    t1 = -2.0 * np.einsum('...ij,...j->...i', At, dalpha)

    # G^i_jk Atilde^jk
    ta = np.einsum('...ijk,...jk->...i', G, At)

    # -(2/3) gtilde^ij d_j K
    gtinv = np.linalg.inv(gt)
    tb = -(2.0/3.0) * np.einsum('...ij,...j->...i', gtinv, dK)

    # 6 Atilde^ij d_j phi
    tc = 6.0 * np.einsum('...ij,...j->...i', At, dphi)

    t2 = 2.0 * alpha[..., None] * (ta + tb + tc)
    return t1 + t2


__all__ = ["dGamma_RHS"]
