"""Sommerfeld outflow BC for all BSSN fields (v2.0-beta3)."""
from __future__ import annotations
import numpy as np


def _apply_scalar(f, f_old, h, dt, c):
    alpha = (c*dt - h) / (c*dt + h)
    out = f.copy()
    for ax in range(3):
        sl = [slice(None)]*3
        sl[ax] = -1
        out[tuple(sl)] = f_old[tuple(sl)] - alpha*(f[tuple(sl)] - f_old[tuple(sl)])
        sl[ax] = 0
        out[tuple(sl)] = f_old[tuple(sl)] - alpha*(f[tuple(sl)] - f_old[tuple(sl)])
    return out


def sommerfeld_all(state, state_old, h, dt, c=1.0):
    """Apply Sommerfeld to phi, K, alpha and each component of gtilde,
    Atilde, Gamma, beta."""
    out = {}
    for key in ("phi", "K", "alpha"):
        out[key] = _apply_scalar(getattr(state, key),
                                  getattr(state_old, key), h, dt, c)
    for key in ("gtilde", "Atilde"):
        arr = getattr(state, key)
        old = getattr(state_old, key)
        new = arr.copy()
        for i in range(3):
            for j in range(3):
                new[..., i, j] = _apply_scalar(arr[..., i, j], old[..., i, j],
                                                h, dt, c)
        out[key] = new
    for key in ("Gamma", "beta"):
        arr = getattr(state, key)
        old = getattr(state_old, key)
        new = arr.copy()
        for i in range(3):
            new[..., i] = _apply_scalar(arr[..., i], old[..., i], h, dt, c)
        out[key] = new
    return out


__all__ = ["sommerfeld_all"]
