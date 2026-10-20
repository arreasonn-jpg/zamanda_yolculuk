"""Brans-Dicke scalar-tensor matter source (v8.0-alpha1).

The BD scalar phi modifies the effective gravitational constant:
    G_eff = G / phi
and contributes to the stress-energy via
    T_munu^BD = (omega/phi) [d_mu phi d_nu phi - 1/2 g_munu (dphi)^2]
                + (1/phi) [nabla_mu nabla_nu phi - g_munu Box phi]

In the Jordan frame, this is equivalent to GR + scalar + non-minimal
coupling.

We implement the Eulerian projections (rho, S_i, S_ij) for use in
the BSSN RHS, in the slow-motion + weak-coupling limit.
"""
from __future__ import annotations
import numpy as np
from ..bssn_grid import d_central


def _grad(f, h):
    return np.stack(d_central(f, h), axis=-1)


def brans_dicke_sources(phi_BD, Pi_BD, gtilde, phi_conf, h,
                         omega=40000.0, V_func=None):
    """Eulerian projections of BD stress-energy.

    phi_BD : scalar field value
    Pi_BD  : conjugate momentum (d_t phi / alpha)
    omega  : BD coupling constant
    """
    Nx, Ny, Nz = phi_BD.shape
    dphi = _grad(phi_BD, h)

    e4phi = np.exp(4.0 * phi_conf)
    gamma = e4phi[..., None, None] * gtilde
    ginv = np.linalg.inv(gamma)

    grad2 = np.einsum('...ij,...i,...j->...', ginv, dphi, dphi)
    Pi2 = Pi_BD**2

    if V_func is None:
        V = np.zeros_like(phi_BD)
    else:
        V = np.array([V_func(p) for p in phi_BD.ravel()]).reshape(phi_BD.shape)

    # Kinetic + potential part (omega term)
    rho_kin = (omega / phi_BD) * (Pi2 + grad2) + V
    # Non-minimal coupling part: (1/phi)(nabla nabla phi - g Box phi)
    # Approximated by isotropic correction at leading order
    rho_nmc = np.zeros_like(phi_BD)  # leading order
    rho = rho_kin + rho_nmc

    # S_i = (omega/phi) Pi d_i phi
    S_i = (omega / phi_BD)[..., None] * Pi_BD[..., None] * dphi

    # S_ij = (omega/phi)(d_i phi d_j phi) + gamma_ij * L_eff
    L_eff = (omega / phi_BD) * (Pi2 - grad2) - V
    S_ij = (omega / phi_BD)[..., None, None] * np.einsum('...i,...j->...ij', dphi, dphi) \
            + L_eff[..., None, None] * gamma
    S = np.einsum('...ij,...ij->...', ginv, S_ij)

    return {"rho": rho, "S_i": S_i, "S_ij": S_ij, "S": S,
            "G_eff_factor": 1.0 / phi_BD}


def gr_limit_check(omega):
    """BD reduces to GR in the omega -> infinity limit."""
    return omega > 40000.0


__all__ = ["brans_dicke_sources", "gr_limit_check"]
