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


def _restrict_into_scalar(f_f, out_c, N_c, offset):
    """Restrict scalar fine -> coarse, overwriting only fully-covered cells."""
    r = restrict(f_f, N_c, offset)
    valid = ~np.isnan(r)
    out_c[valid] = r[valid]
    return out_c


def _restrict_into_tensor(f_f, out_c, N_c, offset):
    """Restrict tensor component fine -> coarse."""
    r = restrict(f_f, N_c, offset)
    valid = ~np.isnan(r)
    out_c[valid] = r[valid]
    return out_c



def restrict_state_into(state_f, coarse_dict, N_c, offset):
    """Restrict fine BSSN arrays into coarse dict (fully-covered cells only)."""
    _restrict_into_scalar(state_f.phi, coarse_dict["phi"], N_c, offset)
    _restrict_into_scalar(state_f.K, coarse_dict["K"], N_c, offset)
    _restrict_into_scalar(state_f.alpha, coarse_dict["alpha"], N_c, offset)
    for i in range(3):
        for j in range(3):
            _restrict_into_tensor(state_f.gtilde[..., i, j],
                                   coarse_dict["gtilde"][..., i, j], N_c, offset)
            _restrict_into_tensor(state_f.Atilde[..., i, j],
                                   coarse_dict["Atilde"][..., i, j], N_c, offset)
    for i in range(3):
        _restrict_into_tensor(state_f.Gamma[..., i],
                               coarse_dict["Gamma"][..., i], N_c, offset)
        _restrict_into_tensor(state_f.beta[..., i],
                               coarse_dict["beta"][..., i], N_c, offset)
    return coarse_dict


__all__ = ["prolongate_state", "restrict_state_into"]
