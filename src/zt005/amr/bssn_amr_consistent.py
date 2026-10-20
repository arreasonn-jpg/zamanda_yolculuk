"""Constraint-consistent AMR boundary injection (v12.15).

Root cause identified in v12.9-v12.13: prolongated boundary values do
not satisfy the fine-grid constraint system. The most directly
constraint-coupled inconsistency is Gamma^i: the stored fine Gamma is
prolongated from coarse (derived with h_c stencils), while the fine RHS
recomputes Christoffel symbols from fine gtilde with h_f stencils every
step. The band Gamma and band gtilde are therefore mutually
inconsistent at exactly the layers where the fine grid is forced.

Fix here: rederive Gamma^i = gtilde^jk G^k_ij from the *prolongated
gtilde* using the fine spacing and the SAME central-difference
operator the RHS uses (bssn_grid.d_central), so that the injected band
is Gamma-consistent with the gtilde it must evolve from.

This is a targeted first step of constraint-consistent prolongation;
it does not touch phi/gtilde/K/Atilde (those mismatches are the known
remaining defect).
"""
from __future__ import annotations

import numpy as np

from ..bssn_grid import d_central


def rederive_Gamma(gtilde, h):
    """Gamma^i = gtilde^jk G^k_ij from gtilde FD derivatives.

    Same Christoffel convention and derivative operator as
    bssn_rhs_vec.bssn_rhs_vec, so the result is consistent with the
    equations actually being integrated.
    """
    gt = np.asarray(gtilde, float)
    dgt = np.zeros(gt.shape + (3,))
    for i in range(3):
        for j in range(3):
            dx, dy, dz = d_central(gt[..., i, j], h)
            dgt[..., i, j, 0] = dx
            dgt[..., i, j, 1] = dy
            dgt[..., i, j, 2] = dz
    gtinv = np.linalg.inv(gt)
    G = 0.5 * (
        np.einsum('...kl,...lji->...kij', gtinv, dgt)
        + np.einsum('...kl,...lij->...kij', gtinv, dgt)
        - np.einsum('...kl,...ijl->...kij', gtinv, dgt))
    # Gamma^i = gtilde^{jk} Gamma^i_{jk}  (upper index FIRST in G storage)
    return np.einsum('...jk,...ijk->...i', gtinv, G)


def gamma_consistent_band(s_f, s_fresh_injected, n_layers, h_f):
    """Replace s_f Gamma in the injected band with FD-derived Gamma.

    s_f / s_fresh_injected : BSSNStateMP-like objects (Gamma attribute).
    Returns Gamma (...,3) with band cells taken from
    rederive_Gamma(s_fresh_injected.gtilde, h_f) and interior kept
    from s_f.Gamma.
    """
    Gamma = np.array(s_f.Gamma, copy=True)
    Gam_f = rederive_Gamma(s_fresh_injected.gtilde, h_f)
    N = Gamma.shape[0]
    for k in range(n_layers):
        Gamma[k, :, :] = Gam_f[k, :, :]
        Gamma[N - 1 - k, :, :] = Gam_f[N - 1 - k, :, :]
        Gamma[:, k, :] = Gam_f[:, k, :]
        Gamma[:, N - 1 - k, :] = Gam_f[:, N - 1 - k, :]
        Gamma[:, :, k] = Gam_f[:, :, k]
        Gamma[:, :, N - 1 - k] = Gam_f[:, :, N - 1 - k]
    return Gamma


__all__ = ["rederive_Gamma", "gamma_consistent_band"]
