"""Puncture-regularized Schwarzschild IC (v2.0-beta6).

Instead of clamping r, use r_p = sqrt(r^2 + eps^2):
  psi(r) = 1 + M / (2 r_p)
  alpha(r) = (1 - M/(2 r_p)) / (1 + M/(2 r_p))

Smooth everywhere. Matches Schwarzschild for r >> eps.
eps = max(M, dx) by default (adaptive).
"""
from __future__ import annotations
import numpy as np
from .bssn_state import BSSNState


def schwarzschild_puncture_state(N, L, M=1.0, eps=None):
    x = np.linspace(-L/2, L/2, N)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    R2 = X**2 + Y**2 + Z**2
    dx = L / (N - 1)
    if eps is None:
        eps = max(M, dx)
    Rp = np.sqrt(R2 + eps**2)
    psi = 1.0 + M / (2.0 * Rp)
    phi = np.log(psi)
    alpha = (1.0 - M/(2.0*Rp)) / (1.0 + M/(2.0*Rp))
    return BSSNState(
        phi=phi,
        gtilde=np.broadcast_to(np.eye(3)[None,None,None,:,:], (N,N,N,3,3)).copy(),
        K=np.zeros((N,N,N)),
        Atilde=np.zeros((N,N,N,3,3)),
        Gamma=np.zeros((N,N,N,3)),
        alpha=alpha,
        beta=np.zeros((N,N,N,3)),
    )


__all__ = ["schwarzschild_puncture_state"]
