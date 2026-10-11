"""Two-level static refinement (AMR-lite, v6.0-alpha3).

Level 0 (coarse): N0 grid, spacing h0, domain [-L/2, L/2]^3
Level 1 (fine):   N1 grid, spacing h1 = h0/2, centered on origin,
                  domain [-L1/2, L1/2]^3 with L1 = L0/2 (typically)

Provides:
- restriction: fine -> coarse (average of 8 fine cells)
- prolongation: coarse -> fine (trilinear interpolation)
- interior index maps
"""
from __future__ import annotations
import numpy as np


def prolongate(f_coarse, N_fine, offset):
    """Trilinear interpolation from coarse grid to fine grid.

    Args:
        f_coarse: (N0, N0, N0)
        N_fine:   int, fine grid size
        offset:   (ox, oy, oz) integer offset: coarse index of fine-cell (0,0,0)
                  in coarse-grid index space.

    Returns:
        f_fine: (N_fine, N_fine, N_fine)
    """
    N0 = f_coarse.shape[0]
    ox, oy, oz = offset
    # Fine cell i corresponds to coarse position ox + i/2
    xs = ox + np.arange(N_fine) / 2.0
    # Clip to valid interior
    xs = np.clip(xs, 0, N0 - 1)
    i0 = np.floor(xs).astype(int)
    i1 = np.minimum(i0 + 1, N0 - 1)
    w = xs - i0
    # Apply along each axis
    def interp_axis(f, idx0, idx1, w, axis):
        f0 = np.take(f, idx0, axis=axis)
        f1 = np.take(f, idx1, axis=axis)
        shape = [1, 1, 1]
        shape[axis] = len(w)
        return f0 * (1 - w.reshape(shape)) + f1 * w.reshape(shape)
    g = interp_axis(f_coarse, i0, i1, w, 0)
    g = interp_axis(g, i0, i1, w, 1)
    g = interp_axis(g, i0, i1, w, 2)
    return g


def restrict(f_fine, N_coarse, offset):
    """Simple averaging restriction: fine -> coarse.

    For each coarse cell (i,j,k), average the 8 fine cells that
    correspond to it. Cells without fine coverage are filled with
    the nearest fine value (fallback).
    """
    N_f = f_fine.shape[0]
    ox, oy, oz = offset
    out = np.zeros((N_coarse, N_coarse, N_coarse))
    for i in range(N_coarse):
        for j in range(N_coarse):
            for k in range(N_coarse):
                # Fine indices covering coarse cell (i,j,k):
                fi0 = 2*(i - ox); fi1 = fi0 + 1
                fj0 = 2*(j - oy); fj1 = fj0 + 1
                fk0 = 2*(k - oz); fk1 = fk0 + 1
                # Fallback: nearest fine
                fii = min(max(fi0, 0), N_f - 1)
                fjj = min(max(fj0, 0), N_f - 1)
                fkk = min(max(fk0, 0), N_f - 1)
                if (0 <= fi0 < N_f - 1 and 0 <= fj0 < N_f - 1
                        and 0 <= fk0 < N_f - 1):
                    patch = f_fine[fi0:fi1+1, fj0:fj1+1, fk0:fk1+1]
                    out[i,j,k] = patch.mean()
                else:
                    out[i,j,k] = f_fine[fii, fjj, fkk]
    return out


__all__ = ["prolongate", "restrict"]
