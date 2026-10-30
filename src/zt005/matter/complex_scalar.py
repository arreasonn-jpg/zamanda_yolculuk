"""Complex scalar matter module (v6.0-alpha1).

Field content:
  psi = psi_R + i psi_I  : complex scalar field
  Pi_psi = Pi_R + i Pi_I : conjugate momentum

U(1) symmetric potential:
  V(psi) = m^2 |psi|^2 + lambda |psi|^4

Stress-energy tensor:
  T_mu_nu = D_mu psi* D_nu psi + D_nu psi* D_mu psi
            - g_mu_nu (D^alpha psi* D_alpha psi + V)

Eulerian projections (3+1 split):
  rho = n^mu n^nu T_mu_nu
  S_i = -gamma_i^mu n^nu T_mu_nu
  S_ij = gamma_i^mu gamma_j^nu T_mu_nu
"""
from __future__ import annotations
import numpy as np
from ..bssn_grid import d_central


def _grad(f, h):
    return np.stack(d_central(f, h), axis=-1)


def complex_scalar_sources(psi_R, psi_I, Pi_R, Pi_I, gtilde, phi_conf, h,
                           m2=1.0, lam=1.0):
    """Total rho, S_i, S_ij from complex scalar field.

    psi = psi_R + i psi_I
    Pi = Pi_R + i Pi_I

    Returns dict with rho, S_i, S_ij, S (trace).
    """
    Nx, Ny, Nz = psi_R.shape

    dpsi_R = _grad(psi_R, h)
    dpsi_I = _grad(psi_I, h)

    e4phi = np.exp(4.0 * phi_conf)
    gamma = e4phi[..., None, None] * gtilde
    ginv = np.linalg.inv(gamma)

    # |grad psi|^2 = grad psi_R^2 + grad psi_I^2
    grad2_R = np.einsum('...ij,...i,...j->...', ginv, dpsi_R, dpsi_R)
    grad2_I = np.einsum('...ij,...i,...j->...', ginv, dpsi_I, dpsi_I)
    grad2 = grad2_R + grad2_I

    # |psi|^2
    psi2 = psi_R**2 + psi_I**2

    # |Pi|^2
    Pi2 = Pi_R**2 + Pi_I**2

    # Potential V = m^2 |psi|^2 + lambda |psi|^4
    V = m2 * psi2 + lam * psi2**2

    # Energy density
    rho = 0.5 * Pi2 + 0.5 * grad2 + V

    # Lagrangian density
    L_eff = 0.5 * Pi2 - 0.5 * grad2 - V

    # Momentum density S_i = Pi_R grad psi_R + Pi_I grad psi_I
    S_i = Pi_R[..., None] * dpsi_R + Pi_I[..., None] * dpsi_I

    # Stress S_ij = grad_i psi_R grad_j psi_R + grad_i psi_I grad_j psi_I + L_eff gamma_ij
    S_ij = (np.einsum('...i,...j->...ij', dpsi_R, dpsi_R)
            + np.einsum('...i,...j->...ij', dpsi_I, dpsi_I)
            + L_eff[..., None, None] * gamma)

    S = np.einsum('...ij,...ij->...', ginv, S_ij)

    return {"rho": rho, "S_i": S_i, "S_ij": S_ij, "S": S}


def complex_scalar_field_rhs(psi_R, psi_I, Pi_R, Pi_I, alpha, K, h,
                             m2=1.0, lam=1.0):
    """Field equations for complex scalar:
      d_t psi_R = alpha Pi_R
      d_t psi_I = alpha Pi_I
      d_t Pi_R = alpha K Pi_R + alpha Lap(psi_R) - alpha dV/dpsi_R
      d_t Pi_I = alpha K Pi_I + alpha Lap(psi_I) - alpha dV/dpsi_I

    dV/dpsi_R = 2 m^2 psi_R + 4 lambda |psi|^2 psi_R
    dV/dpsi_I = 2 m^2 psi_I + 4 lambda |psi|^2 psi_I
    """
    def lap(f):
        grad = np.stack(d_central(f, h), axis=-1)
        L = np.zeros_like(f)
        for ax in range(3):
            dd = d_central(grad[..., ax], h)
            L += dd[ax]
        return L

    psi2 = psi_R**2 + psi_I**2
    dV_dpsiR = 2.0 * m2 * psi_R + 4.0 * lam * psi2 * psi_R
    dV_dpsiI = 2.0 * m2 * psi_I + 4.0 * lam * psi2 * psi_I

    dpsiR_t = alpha * Pi_R
    dpsiI_t = alpha * Pi_I
    dPiR_t = alpha * K * Pi_R + alpha * lap(psi_R) - alpha * dV_dpsiR
    dPiI_t = alpha * K * Pi_I + alpha * lap(psi_I) - alpha * dV_dpsiI

    return dpsiR_t, dpsiI_t, dPiR_t, dPiI_t


def complex_scalar_profile(N, L, psi0=0.1, sigma=2.0):
    """Gaussian complex scalar profile (real part only, imaginary = 0)."""
    x = np.linspace(-L/2, L/2, N)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    R2 = X**2 + Y**2 + Z**2
    psi_R = psi0 * np.exp(-R2 / (2.0 * sigma**2))
    psi_I = np.zeros_like(psi_R)
    Pi_R = np.zeros_like(psi_R)
    Pi_I = np.zeros_like(psi_R)
    return psi_R, psi_I, Pi_R, Pi_I


__all__ = [
    "complex_scalar_sources",
    "complex_scalar_field_rhs",
    "complex_scalar_profile",
]