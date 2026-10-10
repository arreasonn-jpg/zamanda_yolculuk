"""Slowly-rotating Kerr IC for BSSN (v2.0-beta7).

Strategy: Schwarzschild puncture + frame-dragging shift beta^phi.

For a << M in isotropic-like coords:
  beta^phi ~ 2 a M / r^3  (leading order)
  g_tphi ~ -2 a M sin^2(theta) / r
All other BSSN fields = Schwarzschild puncture.

This is NOT full Boyer-Lindquist Kerr; it is a leading-order
slow-rotation approximation suitable for testing a-dependence.
"""
from __future__ import annotations
import numpy as np
from .bssn_state import BSSNState
from .bssn_ic_puncture import schwarzschild_puncture_state


def kerr_slow_rotation_state(N, L, M=1.0, a=0.1, eps=None):
    if abs(a) > 0.5 * M:
        raise ValueError("a must satisfy |a| <= M/2 for slow-rotation")
    s = schwarzschild_puncture_state(N, L, M=M, eps=eps)
    x = np.linspace(-L/2, L/2, N)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    R2 = X**2 + Y**2 + Z**2
    dx = L / (N - 1)
    if eps is None:
        eps = max(M, dx)
    Rp = np.sqrt(R2 + eps**2)
    # beta^phi ~ 2 a M / r^3 ; convert to Cartesian components
    mag = 2.0 * a * M / (Rp**3)
    # unit vector in phi direction: (-Y, X, 0)/sqrt(X^2+Y^2)
    rho = np.sqrt(X**2 + Y**2)
    rho_safe = np.maximum(rho, 1e-12)
    bx = -Y / rho_safe * mag
    by =  X / rho_safe * mag
    bz = np.zeros_like(bx)
    s.beta[..., 0] = bx
    s.beta[..., 1] = by
    s.beta[..., 2] = bz
    return s


__all__ = ["kerr_slow_rotation_state"]
