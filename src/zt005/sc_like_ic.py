"""SC-like smoothed-core initial data for BSSN (v3.0-alpha5).

A conformally-flat regular core approximation of Esa (2026) Schwarzschild
Condensation:
  M_eff(r) = M (1 - e^{-r/rc} (1 + r/rc))
  psi(r)   = 1 + M_eff(r) / (2 r)
  alpha(r) = (1 - M_eff/(2r)) / (1 + M_eff/(2r))
  gtilde   = I, K = 0, Atilde = 0, Gamma = 0, beta = 0

Properties:
- psi(0) = 1      (no central singularity)
- alpha(0) = 1    (regular lapse)
- r >> rc: M_eff -> M, recovers Schwarzschild isotropic
"""
from __future__ import annotations
import numpy as np
from .bssn_state import BSSNState


def sc_like_state(N, L, M=1.0, rc=1.0, eps=1e-6):
    """Smoothed-core BSSN state."""
    x = np.linspace(-L/2, L/2, N)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    R = np.sqrt(X**2 + Y**2 + Z**2)
    R = np.maximum(R, eps)
    exp_term = np.exp(-R/rc) * (1.0 + R/rc)
    M_eff = M * (1.0 - exp_term)
    psi = 1.0 + M_eff / (2.0 * R)
    phi = np.log(psi)
    num = 1.0 - M_eff / (2.0 * R)
    den = 1.0 + M_eff / (2.0 * R)
    alpha = np.where(den > 1e-10, num / den, 1.0)
    alpha = np.maximum(alpha, 1e-6)
    return BSSNState(
        phi=phi,
        gtilde=np.broadcast_to(np.eye(3)[None,None,None,:,:],
                                (N,N,N,3,3)).copy(),
        K=np.zeros((N,N,N)),
        Atilde=np.zeros((N,N,N,3,3)),
        Gamma=np.zeros((N,N,N,3)),
        alpha=alpha,
        beta=np.zeros((N,N,N,3)),
    )


def mass_function_profile(N, L, M=1.0, rc=1.0):
    x = np.linspace(-L/2, L/2, N)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    R = np.sqrt(X**2 + Y**2 + Z**2) + 1e-6
    M_eff = M * (1.0 - np.exp(-R/rc) * (1.0 + R/rc))
    return R, M_eff


__all__ = ["sc_like_state", "mass_function_profile"]
