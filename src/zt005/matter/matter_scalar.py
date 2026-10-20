"""Real scalar field 3+1 stress-energy sources for BSSN (v4.0-alpha1).

Units: G = c = 1 (geometrized).

Given BSSN variables (alpha, beta, gtilde, phi_conformal) and a scalar
field (phi_field, Pi_field) with Pi = (1/alpha)(d_t phi - beta^i d_i phi),
this module computes the Eulerian projections:

  rho  = 1/2 Pi^2 + 1/2 |grad phi|^2 + V(phi)
  S_i  = Pi * d_i phi
  S_ij = d_i phi d_j phi + gamma_ij * L
  S    = gamma^ij S_ij

where
  L       = 1/2 Pi^2 - 1/2 |grad phi|^2 - V(phi)
  gamma_ij = e^{4 phi_c} gtilde_ij
  |grad phi|^2 = gamma^ij d_i phi d_j phi
"""
from __future__ import annotations
import numpy as np
from ..bssn_grid import d_central


def _grad(f, h):
    return np.stack(d_central(f, h), axis=-1)


def scalar_sources(phi_field, Pi_field, gamma_tilde, phi_conf, h,
                   V_func=None):
    """Return dict with rho, S_i, S_ij, S arrays."""
    Nx, Ny, Nz = phi_field.shape
    dphi = _grad(phi_field, h)            # (N,N,N,3)
    # gamma_ij = e^{4 phi_c} gtilde_ij
    e4phi = np.exp(4.0 * phi_conf)
    gamma = e4phi[..., None, None] * gamma_tilde
    ginv = np.linalg.inv(gamma)

    # |grad phi|^2 = gamma^ij d_i phi d_j phi
    grad2 = np.einsum('...ij,...i,...j->...', ginv, dphi, dphi)

    if V_func is None:
        V = np.zeros_like(phi_field)
    else:
        V = np.array([V_func(p) for p in phi_field.ravel()]).reshape(phi_field.shape)

    rho = 0.5 * Pi_field**2 + 0.5 * grad2 + V
    S_i = Pi_field[..., None] * dphi                  # (N,N,N,3)
    L_eff = 0.5 * Pi_field**2 - 0.5 * grad2 - V
    # S_ij = d_i phi d_j phi + gamma_ij * L
    S_ij = np.einsum('...i,...j->...ij', dphi, dphi) + L_eff[..., None, None] * gamma
    S = np.einsum('...ij,...ij->...', ginv, S_ij)

    return {"rho": rho, "S_i": S_i, "S_ij": S_ij, "S": S}


def hamiltonian_matter_term(rho, alpha):
    """4 pi alpha rho -> term added to dK."""
    return 4.0 * np.pi * alpha * rho


def momentum_matter_term(S_i, alpha, gamma_tilde, phi_conf):
    """-16 pi alpha gtilde^ij S_j."""
    e4phi = np.exp(4.0 * phi_conf)
    ginv = np.linalg.inv(e4phi[..., None, None] * gamma_tilde)
    return -16.0 * np.pi * alpha[..., None] * np.einsum('...ij,...j->...i', ginv, S_i)


def atilde_matter_term(S_ij, S, gamma, alpha, phi_conf):
    """8 pi alpha e^{-4 phi} [S_ij - (1/3) gamma_ij S]^TF."""
    trace_term = (1.0/3.0) * S[..., None, None] * gamma
    STF = S_ij - trace_term
    return 8.0 * np.pi * alpha[..., None, None] * np.exp(-4.0*phi_conf)[..., None, None] * STF


__all__ = ["scalar_sources", "hamiltonian_matter_term",
           "momentum_matter_term", "atilde_matter_term"]
