"""Rotating SC smoothed-core IC (v3.0-alpha6).

Combines the SC smoothed core (sc_like_ic) with frame-dragging shift:

  M_eff(r) = M (1 - e^{-r/rc} (1 + r/rc))
  psi(r)   = 1 + M_eff(r) / (2r)
  alpha(r) = (1 - M_eff/(2r)) / (1 + M_eff/(2r))
  beta^phi ~ 2 a M_eff / r^3   (leading-order slow rotation)

Only valid for |a| <= 0.5 M (slow-rotation regime).
"""
from __future__ import annotations
import numpy as np
from .bssn_state import BSSNState
from .sc_like_ic import sc_like_state


def sc_like_kerr_state(N, L, M=1.0, rc=1.0, a=0.1, eps=1e-6):
    if abs(a) > 0.5 * M:
        raise ValueError("a must satisfy |a| <= M/2 for slow rotation")
    s = sc_like_state(N, L, M=M, rc=rc, eps=eps)
    x = np.linspace(-L/2, L/2, N)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    R = np.sqrt(X**2 + Y**2 + Z**2) + eps
    exp_term = np.exp(-R/rc) * (1.0 + R/rc)
    M_eff = M * (1.0 - exp_term)
    mag = 2.0 * a * M_eff / (R**3)
    rho = np.sqrt(X**2 + Y**2) + eps
    bx = -Y / rho * mag
    by =  X / rho * mag
    bz = np.zeros_like(bx)
    s.beta[..., 0] = bx
    s.beta[..., 1] = by
    s.beta[..., 2] = bz
    return s


__all__ = ["sc_like_kerr_state"]
