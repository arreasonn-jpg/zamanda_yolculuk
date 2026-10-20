"""BSSN RHS with scalar field matter source terms (v4.0-alpha1).

Wraps the existing bssn_rhs_beta1 and adds:
  dK     += 4 pi alpha (rho + S)
  dAt    += 8 pi alpha e^{-4 phi_c} S_ij^TF
  dGamma += -16 pi alpha gtilde^ij S_j
"""
from __future__ import annotations
import numpy as np
from .bssn_full import bssn_rhs_beta1
from .matter.matter_scalar import scalar_sources
from .bssn_grid import d_central


def bssn_rhs_matter(phi, gt, K, At, alpha, beta, phi_field, Pi_field,
                    h, V_func=None):
    """Full RHS with matter sources."""
    rhs = bssn_rhs_beta1(phi, gt, K, At, alpha, h)

    # matter sources (V approx = scalar only for now)
    src = scalar_sources(phi_field, Pi_field, gt, phi, h, V_func)
    rho, S_i, S_ij, S = src["rho"], src["S_i"], src["S_ij"], src["S"]

    e4phi = np.exp(4.0 * phi)
    gamma = e4phi[..., None, None] * gt
    ginv = np.linalg.inv(gamma)

    # dK: +4 pi alpha (rho + S)
    rhs["dK"] = rhs["dK"] + 4.0 * np.pi * alpha * (rho + S)

    # dAt: +8 pi alpha e^{-4 phi} S_ij^TF (in conformal frame)
    S_trace = (1.0/3.0) * S[..., None, None] * gamma
    STF = S_ij - S_trace
    rhs["dAt"] = rhs["dAt"] + 8.0 * np.pi * alpha[..., None, None] * \
                 np.exp(-4.0*phi)[..., None, None] * STF

    # dGamma: -16 pi alpha gtilde^ij S_j
    rhs["dGamma"] = rhs["dGamma"] - 16.0 * np.pi * alpha[..., None] * \
                    np.einsum('...ij,...j->...i', ginv, S_i)

    return rhs, {"rho": rho, "S_i": S_i, "S_ij": S_ij, "S": S}


__all__ = ["bssn_rhs_matter"]
