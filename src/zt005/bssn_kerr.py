"""Kerr metric in Boyer-Lindquist; a=0 limit = Schwarzschild.
We only check tensor/constraint sanity locally, not full BSSN evolution."""
from __future__ import annotations
import numpy as np


def kerr_bl(M: float, a: float, r: float, theta: float):
    """Return (alpha, gamma_ij) in BL coords at (r, theta)."""
    if abs(a) > M:
        raise ValueError("a <= M required")
    Sigma = r**2 + a**2 * np.cos(theta)**2
    Delta = r**2 - 2.0 * M * r + a**2
    if Delta <= 0:
        raise ValueError("r outside horizon (Delta <= 0)")
    rho2 = Sigma
    # Lapse
    A = (r**2 + a**2)**2 - a**2 * Delta * np.sin(theta)**2
    alpha2 = rho2 * Delta / A
    alpha = np.sqrt(max(alpha2, 1e-300))
    # 3-metric (r, theta, phi) diagonal approx for small a
    g_rr = rho2 / Delta
    g_tt = rho2
    g_pp = A * np.sin(theta)**2 / rho2
    gamma = np.diag([g_rr, g_tt, g_pp])
    return alpha, gamma


def kerr_horizon_radius(M: float, a: float) -> float:
    return M + np.sqrt(max(M**2 - a**2, 0.0))


__all__ = ["kerr_bl", "kerr_horizon_radius"]
