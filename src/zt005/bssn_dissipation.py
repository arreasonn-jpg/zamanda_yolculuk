"""Kreiss-Oliger dissipation + simple Sommerfeld outflow boundary."""
from __future__ import annotations
import numpy as np


def ko_dissipation(f, h, eps=0.1):
    """6th-order Kreiss-Oliger operator on a uniform grid.
    Uses standard (-1)^k D_+ D_- stencil; boundaries left untouched."""
    out = np.zeros_like(f)
    for ax in range(3):
        def shift(arr, n):
            return np.roll(arr, n, axis=ax)
        d = (shift(f, -3) - 6*shift(f, -2) + 15*shift(f, -1)
             - 20*f + 15*shift(f, 1) - 6*shift(f, 2) + shift(f, 3))
        # 6th-order: coeff = eps * h^5 / 64 * (D_+ D_-)^3
        out += (-1) ** 3 * eps * h ** 5 / 64.0 * d
    return out


def sommerfeld_boundary(f, f_old, h, dt, c=1.0):
    """Simple 1st-order outflow on outer faces along each axis.
    f_new = f_old - (c dt - h)/(c dt + h) * (f - f_old) at last index."""
    out = f.copy()
    alpha = (c * dt - h) / (c * dt + h)
    for ax in range(3):
        sl = [slice(None)] * 3
        sl[ax] = -1
        out[tuple(sl)] = f_old[tuple(sl)] - alpha * (f[tuple(sl)] - f_old[tuple(sl)])
    return out


def apply_dissipation_to_state(state, h, dt, eps=0.1):
    """Apply KO dissipation to phi, gtilde, K, Atilde with dt scaling."""
    scale = dt / h
    out = dict(state)
    for key in ("phi", "gtilde", "K", "Atilde"):
        out[key] = state[key] + scale * ko_dissipation(state[key], h, eps)
    return out


__all__ = ["ko_dissipation", "sommerfeld_boundary", "apply_dissipation_to_state"]
