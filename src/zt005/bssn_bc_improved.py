"""Improved boundary conditions with constraint damping (v12.7-alpha1).

Standard Sommerfeld BC only handles the outgoing wave; it does not
suppress constraint violation that accumulates from wave reflection.
We add two improvements:

1. Constraint damping term near boundaries:
       dK += -kappa * H * chi_boundary(x)
   where chi_boundary is a smooth taper (0 in interior, 1 at edge).
   This is a Z4c-inspired sponge layer.

2. Multi-layer Sommerfeld: apply the Sommerfeld update over the
   last N_boundary cells (not just the outermost), with an
   exponential weight.
"""
from __future__ import annotations
import numpy as np

def sponge_profile(N, width_cells):
    """Smooth sponge taper in the outermost `width_cells` layers.
    Returns array (N,) with 0 in interior, 1 at edges.
    """
    prof = np.zeros(N)
    for i in range(N):
        d = min(i, N - 1 - i)
        if d < width_cells:
            x = (width_cells - d) / width_cells  # 0 at inner edge, 1 at outer
            prof[i] = 3 * x**2 - 2 * x**3  # smoothstep
    return prof

def sponge_3d(N, width_cells):
    """3D sponge profile (product of 1D tapers)."""
    p = sponge_profile(N, width_cells)
    Px, Py, Pz = np.meshgrid(p, p, p, indexing="ij")
    # Union: take max (any boundary surface)
    return np.maximum(np.maximum(Px, Py), Pz)

def damp_K(K, H, sponge, kappa=0.5):
    """Apply constraint damping to K near boundaries."""
    return K - kappa * sponge * H

def multi_layer_sommerfeld(f, f_old, h, dt, width_cells, c=1.0):
    """Apply Sommerfeld update over multiple boundary layers."""
    N = f.shape[0]
    out = f.copy()
    # 1D alpha per axis
    alpha = (c * dt - h) / (c * dt + h)
    prof = sponge_profile(N, width_cells)
    for ax in range(3):
        # Build multi-dim sponge by broadcasting
        shape = [1, 1, 1]
        shape[ax] = N
        w = prof.reshape(shape)
        # d_t f + c d_n f = 0 -> f_new = f_old - alpha * (f - f_old) near boundary
        out = out - w * alpha * (f - f_old)
    return out

__all__ = ["sponge_profile", "sponge_3d", "damp_K",
           "multi_layer_sommerfeld"]
