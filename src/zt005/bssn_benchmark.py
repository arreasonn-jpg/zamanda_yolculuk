"""Known-metric benchmarks for the BSSN skeleton.

Schwarzschild in isotropic coords:
  ds^2 = -alpha^2 dt^2 + psi^4 (dr^2 + r^2 dOmega^2)
  alpha = (1 - M/(2r)) / (1 + M/(2r))
  psi   = 1 + M/(2r)
So: gtilde_ij = delta_ij  (conformally flat)
    phi = ln(psi)
    K = 0, Atilde = 0
Vacuum Hamiltonian: H = R + K^2 - Kij K^ij - 16 pi G rho / c^2 = R
For Schwarzschild: R = 0 outside horizon.
"""
from __future__ import annotations
import numpy as np

C = 2.99792458e8
G = 6.67430e-11


def schwarzschild_isotropic(M: float, r: float):
    """Return BSSN (alpha, phi, gtilde, K, Atilde) at radius r."""
    if r <= M / 2.0:
        raise ValueError("r must be > M/2 (isotropic horizon)")
    alpha = (1.0 - M / (2.0 * r)) / (1.0 + M / (2.0 * r))
    psi = 1.0 + M / (2.0 * r)
    return {
        "alpha": alpha,
        "phi": float(np.log(psi)),
        "gtilde": np.eye(3),
        "K": 0.0,
        "Atilde": np.zeros((3, 3)),
    }


def _d_central_1d(f, h):
    d = np.zeros_like(f)
    d[1:-1] = (f[2:] - f[:-2]) / (2.0 * h)
    return d


def ricci_scalar_radial(phi_r, r, h):
    """3D Ricci scalar of gamma_ij = e^{4 phi} delta_ij on a radial line,
    spherical symmetry assumed. Uses FD on phi(r).

    R = -8 e^{-4 phi} (phi'' + (2/r) phi' + (phi')^2)  ... contracted
    Actually for conformally flat 3-metric with psi=e^phi:
      R = -8 e^{-4 phi} (nabla^2 phi + (1/2) |grad phi|^2)
    We compute on the radial line assuming phi=phi(r)."""
    phi_r = np.asarray(phi_r, float)
    phi_p = _d_central_1d(phi_r, h)
    phi_pp = _d_central_1d(phi_p, h)
    lap = phi_pp + 2.0 * phi_p / np.maximum(r, 1e-12)
    grad2 = phi_p ** 2
    return -8.0 * np.exp(-4.0 * phi_r) * (lap + grad2)
