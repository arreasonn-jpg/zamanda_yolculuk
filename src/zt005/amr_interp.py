"""AMR interpolation operators (v12.9-alpha1).

Trilinear prolongation (coarse -> fine) and restriction (fine -> coarse).
Pure geometry, no BSSN physics. Used by fixed-patch AMR in v12.9-alpha2.
"""
from __future__ import annotations
import numpy as np


def prolongate_scalar(coarse, factor=2):
    """Coarse (N,N,N) -> fine. Trilinear interpolation."""
    N = coarse.shape[0]
    Nf = (N - 1) * factor + 1
    xf = np.linspace(0, N - 1, Nf)
    i0 = np.floor(xf).astype(int)
    i1 = np.minimum(i0 + 1, N - 1)
    t = xf - i0
    out = coarse.copy()
    for axis in range(3):
        a = np.take(out, i0, axis=axis)
        b = np.take(out, i1, axis=axis)
        sh = [1, 1, 1]
        sh[axis] = -1
        t_bc = t.reshape(sh)
        out = a * (1 - t_bc) + b * t_bc
    return out


def prolongate_tensor(coarse, factor=2):
    """Coarse (N,N,N,3,3) -> fine. Apply scalar prolongation per component."""
    N = coarse.shape[0]
    Nf = (N - 1) * factor + 1
    out = np.zeros((Nf, Nf, Nf, 3, 3))
    for i in range(3):
        for j in range(3):
            out[..., i, j] = prolongate_scalar(coarse[..., i, j], factor)
    return out


def prolongate_vector(coarse, factor=2):
    """Coarse (N,N,N,3) -> fine."""
    N = coarse.shape[0]
    Nf = (N - 1) * factor + 1
    out = np.zeros((Nf, Nf, Nf, 3))
    for i in range(3):
        out[..., i] = prolongate_scalar(coarse[..., i], factor)
    return out


def restrict_scalar(fine, factor=2):
    """Fine -> coarse by sampling every 'factor' cells (injection).

    Simple injection (not averaging) - adequate for fixed-patch AMR.
    """
    return fine[::factor, ::factor, ::factor].copy()


def inject_boundary(fine, coarse_prolongated, n_layers=4):
    """Overwrite fine boundary with prolongated coarse values.

    Keeps fine interior, replaces n_layers outermost cells with coarse BC.
    """
    Nf = fine.shape[0]
    out = fine.copy()
    # Copy boundary from coarse_prolongated
    for k in range(n_layers):
        # Faces: 6 sides
        out[k, :, :] = coarse_prolongated[k, :, :]
        out[-1-k, :, :] = coarse_prolongated[-1-k, :, :]
        out[:, k, :] = coarse_prolongated[:, k, :]
        out[:, -1-k, :] = coarse_prolongated[:, -1-k, :]
        out[:, :, k] = coarse_prolongated[:, :, k]
        out[:, :, -1-k] = coarse_prolongated[:, :, -1-k]
    return out


__all__ = ["prolongate_scalar", "prolongate_tensor", "prolongate_vector",
           "restrict_scalar", "inject_boundary"]
