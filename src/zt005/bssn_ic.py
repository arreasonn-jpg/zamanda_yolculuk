"""Constraint-satisfying initial data for BSSN (v2.0-beta5).

Schwarzschild in isotropic coordinates (maximally sliced):
  phi(r) = ln(1 + M/(2r))
  K = 0
  Atilde = 0
  gtilde = I
  Gamma = 0
  alpha(r) = (1 - M/(2r)) / (1 + M/(2r))
  beta = 0
"""
from __future__ import annotations
import numpy as np
from .bssn_state import BSSNState


def schwarzschild_isotropic_state(N, L, M=1.0):
    x = np.linspace(-L/2, L/2, N)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    R = np.sqrt(X**2 + Y**2 + Z**2)
    R = np.maximum(R, 0.51 * M)  # avoid r < M/2
    phi = np.log(1.0 + M / (2.0 * R))
    alpha = (1.0 - M/(2.0*R)) / (1.0 + M/(2.0*R))
    return BSSNState(
        phi=phi,
        gtilde=np.broadcast_to(np.eye(3)[None,None,None,:,:], (N,N,N,3,3)).copy(),
        K=np.zeros((N,N,N)),
        Atilde=np.zeros((N,N,N,3,3)),
        Gamma=np.zeros((N,N,N,3)),
        alpha=alpha,
        beta=np.zeros((N,N,N,3)),
    )


__all__ = ["schwarzschild_isotropic_state"]
