"""Tensor-aware prolongate/restrict for BSSN fields (v6.0-alpha5).

Wraps scalar refinement operators for:
- scalar fields: phi, K, alpha
- vector fields: Gamma, beta
- tensor fields: gtilde, Atilde
"""
from __future__ import annotations
import numpy as np
from .refinement import prolongate, restrict


def prolongate_state(state_c, N_f, offset):
    """Prolongate a coarse BSSNState to a fine grid (returns dict of arrays)."""
    out = {}
    out["phi"] = prolongate(state_c.phi, N_f, offset)
    out["K"] = prolongate(state_c.K, N_f, offset)
    out["alpha"] = prolongate(state_c.alpha, N_f, offset)
    # tensor fields
    g = np.zeros((N_f, N_f, N_f, 3, 3))
    A = np.zeros((N_f, N_f, N_f, 3, 3))
    for i in range(3):
        for j in range(3):
            g[..., i, j] = prolongate(state_c.gtilde[..., i, j], N_f, offset)
            A[..., i, j] = prolongate(state_c.Atilde[..., i, j], N_f, offset)
    out["gtilde"] = g
    out["Atilde"] = A
    # vector fields
    G = np.zeros((N_f, N_f, N_f, 3))
    b = np.zeros((N_f, N_f, N_f, 3))
    for i in range(3):
        G[..., i] = prolongate(state_c.Gamma[..., i], N_f, offset)
        b[..., i] = prolongate(state_c.beta[..., i], N_f, offset)
    out["Gamma"] = G
    out["beta"] = b
    return out


def restrict_state_into(state_f, coarse_dict, N_c, offset):
    """Restrict fine BSSN arrays into coarse dict in the fully-covered region."""
    ox, oy, oz = offset
    N_f = state_f.phi.shape[0]
    half = N_f // 2
    a, b = ox, ox + half - 1
    sl = (slice(a, b+1), slice(a, b+1), slice(a, b+1))

    coarse_dict["phi"][sl] = restrict(state_f.phi, N_c, offset)[sl]
    coarse_dict["K"][sl] = restrict(state_f.K, N_c, offset)[sl]
    coarse_dict["alpha"][sl] = restrict(state_f.alpha, N_c, offset)[sl]

    # tensor fields: restrict each component
    for i in range(3):
        for j in range(3):
            coarse_dict["gtilde"][..., i, j][sl] = restrict(
                state_f.gtilde[..., i, j], N_c, offset)[sl]
            coarse_dict["Atilde"][..., i, j][sl] = restrict(
                state_f.Atilde[..., i, j], N_c, offset)[sl]
    for i in range(3):
        coarse_dict["Gamma"][..., i][sl] = restrict(
            state_f.Gamma[..., i], N_c, offset)[sl]
        coarse_dict["beta"][..., i][sl] = restrict(
            state_f.beta[..., i], N_c, offset)[sl]
    return coarse_dict


__all__ = ["prolongate_state", "restrict_state_into"]
