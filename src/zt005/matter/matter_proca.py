"""Proca (massive vector) field matter source (v10.2-alpha1).
Lagrangian: L = -1/4 F_mu_nu F^mu_nu - 1/2 m^2 A_mu A^mu
F_mu_nu = d_mu A_nu - d_nu A_mu

Stress-energy:
 T_mu_nu = F_mu_alpha F_nu^alpha + m^2 A_mu A_nu + g_mu_nu L

Key physics: unlike a scalar field, a vector field has
directional (anisotropic) spatial stress T_ij. The TT
projection is generically nonzero, giving real GW emission.
"""
from __future__ import annotations
import numpy as np
from ..bssn_grid import d_central

def _grad(f, h):
    return np.stack(d_central(f, h), axis=-1)

def proca_sources(A0, A1, A2, A3,
                  E1, E2, E3,
                  gtilde, phi_conf, h,
                  m2=0.0):
    """Eulerian projections of Proca stress-energy.

    A0 : scalar potential (A_t)
    A1, A2, A3 : spatial vector components
    E_i : conjugate momenta  (E_i = d_t A_i - d_i A_0)
    m2 : mass squared (0 for Maxwell)
    """
    Nx, Ny, Nz = A0.shape
    # Spatial gradients of A_i (for F_ij = d_i A_j - d_j A_i)
    dA1 = _grad(A1, h)
    dA2 = _grad(A2, h)
    dA3 = _grad(A3, h)

    e4phi = np.exp(4.0 * phi_conf)
    gamma = e4phi[..., None, None] * gtilde
    ginv = np.linalg.inv(gamma)

    # Electric field magnitude squared
    Emag2 = E1**2 + E2**2 + E3**2

    # Magnetic field: B_i = epsilon_ijk d_j A_k  (flat space)
    B1 = dA3[..., 1] - dA2[..., 2]
    B2 = dA1[..., 2] - dA3[..., 0]
    B3 = dA2[..., 0] - dA1[..., 1]
    B2_mag = B1**2 + B2**2 + B3**2

    # Proca energy density
    rho = 0.5 * (Emag2 + B2_mag) + 0.5 * m2 * (A1**2 + A2**2 + A3**2)

    # Poynting vector  S_i = (E x B)_i
    S1 = E2 * B3 - E3 * B2
    S2 = E3 * B1 - E1 * B3
    S3 = E1 * B2 - E2 * B1
    S_i = np.stack([S1, S2, S3], axis=-1)

    # Anisotropic spatial stress (Maxwell + mass term)
    # T_ij = -E_i E_j - B_i B_j + 0.5 delta_ij (E^2 + B^2) + m^2 A_i A_j
    A_vec = np.stack([A1, A2, A3], axis=-1)
    E_vec = np.stack([E1, E2, E3], axis=-1)
    B_vec = np.stack([B1, B2, B3], axis=-1)
    EE = np.einsum('...i,...j->...ij', E_vec, E_vec)
    BB = np.einsum('...i,...j->...ij', B_vec, B_vec)
    AA = np.einsum('...i,...j->...ij', A_vec, A_vec)

    L_eff = 0.5 * (Emag2 - B2_mag) - 0.5 * m2 * (A1**2 + A2**2 + A3**2)
    S_ij = (-EE - BB + m2 * AA
            + L_eff[..., None, None] * gamma)

    S = np.einsum('...ij,...ij->...', ginv, S_ij)
    return {"rho": rho, "S_i": S_i, "S_ij": S_ij, "S": S}

def tt_content_estimate(S_ij, gamma, k_hat, e1, e2):
    """Rough TT projection of S_ij to check GW-relevant content."""
    ginv = np.linalg.inv(gamma)
    # Decompose using projector
    S_mixed = np.einsum('...ij,...jk->...ik', ginv, S_ij)
    # Just return the traceless part norm as proxy
    tr = np.einsum('...ij,...ij->...', ginv, S_ij) / 3.0
    STF = S_ij - tr[..., None, None] * gamma
    return np.sqrt(np.einsum('...ij,...ij->...', STF, STF))

__all__ = ["proca_sources", "tt_content_estimate"]
