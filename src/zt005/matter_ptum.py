"""PTUM matter coupling to BSSN (v5.0-alpha2).

Field content:
  Phi (fuel)    : real scalar
  Pi_Phi        : conjugate momentum
  psi (proxy)   : real scalar
  Pi_psi        : conjugate momentum

Potential:
  V = V_fuel(Phi) + V_proxy(psi; Phi)

Eulerian projections (same structure as matter_scalar.py but for two
fields).
"""
from __future__ import annotations
import numpy as np
from .bssn_grid import d_central
from .ptum_field import V_fuel, V_proxy


def _grad(f, h):
    return np.stack(d_central(f, h), axis=-1)


def ptum_sources(Phi, Pi_Phi, psi, Pi_psi, gtilde, phi_conf, h,
                 V0=1.0, lam=0.05, m2=1.0, alpha=1.0, Eth=1.0,
                 lam_psi=1.0):
    """Total rho, S_i, S_ij from both PTUM fields."""
    Nx, Ny, Nz = Phi.shape
    dPhi = _grad(Phi, h)
    dpsi = _grad(psi, h)

    e4phi = np.exp(4.0 * phi_conf)
    gamma = e4phi[..., None, None] * gtilde
    ginv = np.linalg.inv(gamma)

    grad2_Phi = np.einsum('...ij,...i,...j->...', ginv, dPhi, dPhi)
    grad2_psi = np.einsum('...ij,...i,...j->...', ginv, dpsi, dpsi)

    Vf = V_fuel(Phi, V0=V0, lam=lam)
    Vp = V_proxy(psi, Phi, m2=m2, alpha=alpha, Eth=Eth, lam_psi=lam_psi)
    V_tot = Vf + Vp

    rho = (0.5 * Pi_Phi**2 + 0.5 * Pi_psi**2
           + 0.5 * grad2_Phi + 0.5 * grad2_psi + V_tot)
    L_eff = (0.5 * Pi_Phi**2 + 0.5 * Pi_psi**2
             - 0.5 * grad2_Phi - 0.5 * grad2_psi - V_tot)
    S_i = Pi_Phi[..., None] * dPhi + Pi_psi[..., None] * dpsi
    S_ij = (np.einsum('...i,...j->...ij', dPhi, dPhi)
            + np.einsum('...i,...j->...ij', dpsi, dpsi)
            + L_eff[..., None, None] * gamma)
    S = np.einsum('...ij,...ij->...', ginv, S_ij)

    return {"rho": rho, "S_i": S_i, "S_ij": S_ij, "S": S}


__all__ = ["ptum_sources"]
