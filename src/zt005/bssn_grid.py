"""Finite-difference derivatives on a uniform 3D grid, for BSSN RHS."""
from __future__ import annotations
import numpy as np


def d_central(f, h):
    """Central 2nd-order derivative along each of 3 axes.
    f: (Nx,Ny,Nz). Returns (dfx, dfy, dfz) same shape (zero at edges)."""
    out = []
    for ax in range(3):
        d = np.zeros_like(f)
        sl_c = [slice(None)] * 3
        sl_p = [slice(None)] * 3
        sl_m = [slice(None)] * 3
        sl_c[ax] = slice(1, -1); sl_p[ax] = slice(2, None); sl_m[ax] = slice(None, -2)
        d[tuple(sl_c)] = (f[tuple(sl_p)] - f[tuple(sl_m)]) / (2.0 * h)
        out.append(d)
    return tuple(out)


def dgtilde_from_grid(gtilde, h):
    """gtilde: (Nx,Ny,Nz,3,3) -> dgtilde[...,i,j,k] = d_k gtilde_ij."""
    Nx, Ny, Nz = gtilde.shape[:3]
    dg = np.zeros((Nx, Ny, Nz, 3, 3, 3))
    for i in range(3):
        for j in range(3):
            comp = gtilde[..., i, j]
            dx, dy, dz = d_central(comp, h)
            dg[..., i, j, 0] = dx
            dg[..., i, j, 1] = dy
            dg[..., i, j, 2] = dz
    return dg


__all__ = ["d_central", "dgtilde_from_grid"]
