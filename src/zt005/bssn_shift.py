"""Shift-vector Lie derivatives for BSSN v2.0-beta2."""
from __future__ import annotations
import numpy as np
from .bssn_grid import d_central


def _grad(f, h):
    return np.stack(d_central(f, h), axis=-1)


def _grad3(A, h):
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


def shift_terms(phi, gt, K, At, Gamma, beta, h):
    """Lie-derivative terms along beta for all BSSN variables."""
    dphi = _grad(phi, h)
    dK = _grad(K, h)
    dgt = _grad3(gt, h)
    dAt = _grad3(At, h)
    dbeta = _grad(beta, h)
    dGamma = _grad(Gamma, h)

    Nx, Ny, Nz = phi.shape

    Lphi = np.einsum('...k,...k->...', beta, dphi)
    LK = np.einsum('...k,...k->...', beta, dK)
    Lgt = np.einsum('...k,...ijk->...ij', beta, dgt)
    LAt = np.einsum('...k,...ijk->...ij', beta, dAt)

    divb = dbeta[..., 0, 0] + dbeta[..., 1, 1] + dbeta[..., 2, 2]

    # +2 gtilde_k(i d_j) beta^k
    term_gt = np.zeros_like(gt)
    for i in range(3):
        for j in range(3):
            for k in range(3):
                term_gt[..., i, j] += gt[..., k, i] * dbeta[..., k, j]
                term_gt[..., i, j] += gt[..., k, j] * dbeta[..., k, i]
    Lgt = Lgt + term_gt

    # +2 Atilde_k(i d_j) beta^k - (2/3) Atilde_ij div(beta)
    term_At = np.zeros_like(At)
    for i in range(3):
        for j in range(3):
            for k in range(3):
                term_At[..., i, j] += At[..., k, i] * dbeta[..., k, j]
                term_At[..., i, j] += At[..., k, j] * dbeta[..., k, i]
    LAt = LAt + term_At - (2.0/3.0) * At * divb[..., None, None]

    # Gamma^i shift terms
    LGamma = np.einsum('...j,...ij->...i', beta, dGamma)
    LGamma = LGamma - np.einsum('...j,...ji->...i', Gamma, dbeta)
    divb_grad = np.stack(d_central(divb, h), axis=-1)
    gtinv = np.linalg.inv(gt)
    LGamma = LGamma + (1.0/3.0) * np.einsum('...ij,...j->...i', gtinv, divb_grad)
    LGamma = LGamma - (2.0/3.0) * Gamma * divb[..., None]

    return {"Lphi": Lphi, "Lgt": Lgt, "LK": LK, "LAt": LAt, "LGamma": LGamma}


__all__ = ["shift_terms"]
