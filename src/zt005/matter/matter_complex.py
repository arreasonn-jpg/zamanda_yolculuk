"""Complex scalar field stress-energy for BSSN (v6.0-alpha1).

Complex scalar Psi = Psi_R + i Psi_I (two real fields).

Lagrangian:
  L = -g^mu_nu d_mu Psi* d_nu Psi - V(|Psi|^2)

Eulerian projections:
  rho  = |Pi|^2 + |grad Psi|^2 + V
  S_i  = 2 Re(Pi* d_i Psi)
  S_ij = 2 Re(d_i Psi* d_j Psi) + gamma_ij L_eff
  L_eff = |Pi|^2 - |grad Psi|^2 - V

where
  Pi = (1/alpha)(d_t Psi - beta^i d_i Psi)   (conjugate momentum)
  |Pi|^2 = Pi_R^2 + Pi_I^2
  |grad Psi|^2 = gamma^ij (d_i Psi_R d_j Psi_R + d_i Psi_I d_j Psi_I)
  V = V(|Psi|^2),  |Psi|^2 = Psi_R^2 + Psi_I^2
"""
from __future__ import annotations
import numpy as np
from ..bssn_grid import d_central


def _grad(f, h):
    return np.stack(d_central(f, h), axis=-1)


def complex_sources(Psi_R, Psi_I, Pi_R, Pi_I, gtilde, phi_conf, h,
                    V_func=None):
    """Return dict with rho, S_i, S_ij, S for complex scalar."""
    Nx, Ny, Nz = Psi_R.shape
    dR = _grad(Psi_R, h)
    dI = _grad(Psi_I, h)

    e4phi = np.exp(4.0 * phi_conf)
    gamma = e4phi[..., None, None] * gtilde
    ginv = np.linalg.inv(gamma)

    # |grad Psi|^2
    gradR2 = np.einsum('...ij,...i,...j->...', ginv, dR, dR)
    gradI2 = np.einsum('...ij,...i,...j->...', ginv, dI, dI)
    grad2 = gradR2 + gradI2

    # |Pi|^2
    Pi2 = Pi_R**2 + Pi_I**2

    # |Psi|^2
    Psi2 = Psi_R**2 + Psi_I**2

    if V_func is None:
        V = np.zeros_like(Psi_R)
    else:
        V = np.array([V_func(p) for p in Psi2.ravel()]).reshape(Psi_R.shape)

    rho = Pi2 + grad2 + V
    L_eff = Pi2 - grad2 - V

    # S_i = 2 Re(Pi* d_i Psi) = 2(Pi_R d_i Psi_R + Pi_I d_i Psi_I)
    S_i = 2.0 * (Pi_R[..., None] * dR + Pi_I[..., None] * dI)

    # S_ij = 2 Re(d_i Psi* d_j Psi) + gamma_ij L_eff
    dPsi_dPsi = (np.einsum('...i,...j->...ij', dR, dR)
                 + np.einsum('...i,...j->...ij', dI, dI))
    S_ij = 2.0 * dPsi_dPsi + L_eff[..., None, None] * gamma
    S = np.einsum('...ij,...ij->...', ginv, S_ij)

    return {"rho": rho, "S_i": S_i, "S_ij": S_ij, "S": S}


def charge_density(Psi_R, Psi_I, Pi_R, Pi_I):
    """U(1) charge density: rho_Q = i(Psi* Pi - Pi* Psi) = 2(Psi_R Pi_I - Psi_I Pi_R)."""
    return 2.0 * (Psi_R * Pi_I - Psi_I * Pi_R)


__all__ = ["complex_sources", "charge_density"]
