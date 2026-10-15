"""Fully vectorized BSSN RHS (v12.2-alpha1).

Replaces Python point loops with numpy einsum. Same equations,
~5-10x speedup at N>=11.
"""
from __future__ import annotations
import numpy as np
from .bssn_grid import d_central

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

def _grad_tensor(A, h):
    """A: (...,3,3). Returns dA[...,i,j,k] = d_k A_ij."""
    dA = np.zeros(A.shape + (3,))
    for i in range(3):
        for j in range(3):
            dx, dy, dz = d_central(A[..., i, j], h)
            dA[..., i, j, 0] = dx
            dA[..., i, j, 1] = dy
            dA[..., i, j, 2] = dz
    return dA

def _grad_tensor3(G, h):
    """G: (...,3,3,3). Returns dG[...,k,i,j,l] = d_l Gamma^k_ij."""
    dG = np.zeros(G.shape + (3,))
    for k in range(3):
        for i in range(3):
            for j in range(3):
                dx, dy, dz = d_central(G[..., k, i, j], h)
                dG[..., k, i, j, 0] = dx
                dG[..., k, i, j, 1] = dy
                dG[..., k, i, j, 2] = dz
    return dG

def bssn_rhs_vec(phi, gt, K, At, alpha, h):
    """Vectorized BSSN RHS (alpha=1, beta=0)."""
    dphi = _grad(phi, h)
    Hphi = _hess(phi, h)
    dalpha = _grad(alpha, h)
    Halpha = _hess(alpha, h)
    dK = _grad(K, h)
    dgt = _grad_tensor(gt, h)

    gtinv = np.linalg.inv(gt)

    # Christoffel
    G = 0.5 * (
        np.einsum('...kl,...lji->...kij', gtinv, dgt)
        + np.einsum('...kl,...lij->...kij', gtinv, dgt)
        - np.einsum('...kl,...ijl->...kij', gtinv, dgt)
    )

    # Ricci
    dG = _grad_tensor3(G, h)
    R_ij = (
        np.einsum('...kijk->...ij', dG)
        - np.einsum('...kikj->...ij', dG)
        + np.einsum('...kkl,...lij->...ij', G, G)
        - np.einsum('...kil,...ljk->...ij', G, G)
    )

    # Atilde^2
    A2 = np.einsum('...ik,...kl,...lj->...ij', At, gtinv, At)
    A2_scalar = np.einsum('...ij,...ij->...', At, np.einsum('...ik,...kl->...il', At, gtinv))

    # dK
    DDalpha = (
        np.einsum('...ij,...ij->...', gtinv, Halpha)
        - np.einsum('...ij,...kij,...k->...', gtinv, G, dalpha)
    )
    dK_t = -DDalpha + alpha * (A2_scalar + K**2 / 3.0)

    # dphi
    dphi_t = -alpha * K / 6.0

    # dgt
    dgt_t = -2.0 * alpha[..., None, None] * At

    # dAtilde
    DDphi = Hphi - np.einsum('...kij,...k->...ij', G, dphi)
    dphidphi = np.einsum('...i,...j->...ij', dphi, dphi)
    em4phi = np.exp(-4.0 * phi)
    bracket = (-2.0 * alpha[..., None, None] * DDphi
               + 4.0 * alpha[..., None, None] * dphidphi
               + alpha[..., None, None] * R_ij)
    tr_b = np.einsum('...ij,...ij->...', gtinv, bracket)
    bracket_tf = bracket - (tr_b / 3.0)[..., None, None] * gt
    dAt_t = (em4phi[..., None, None] * bracket_tf
             + alpha[..., None, None] * (K[..., None, None] * At - 2.0 * A2))
    tr_dAt = np.einsum('...ij,...ij->...', gtinv, dAt_t)
    dAt_t = dAt_t - (tr_dAt / 3.0)[..., None, None] * gt

    # dGamma
    At_up = np.einsum('...ik,...kl,...lj->...ij', gtinv, At, gtinv)
    t1 = -2.0 * np.einsum('...ij,...j->...i', At_up, dalpha)
    ta = np.einsum('...ijk,...jk->...i', G, At_up)
    tb = -(2.0/3.0) * np.einsum('...ij,...j->...i', gtinv, dK)
    tc = 6.0 * np.einsum('...ij,...j->...i', At_up, dphi)
    dG_t = t1 + 2.0 * alpha[..., None] * (ta + tb + tc)

    return {"dphi": dphi_t, "dgt": dgt_t, "dK": dK_t,
            "dAt": dAt_t, "dGamma": dG_t}

__all__ = ["bssn_rhs_vec"]
