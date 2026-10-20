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

    For each coarse cell (i,j,k) covered by the fine region, average
    the 2x2x2 fine cells that correspond to it. Cells OUTSIDE the
    fine region are filled with NaN — caller must not use them.
    """
    N_f = f_fine.shape[0]
    ox, oy, oz = offset
    out = np.full((N_coarse, N_coarse, N_coarse), np.nan)
    # Number of coarse cells fully covered: fine point 2*i corresponds to
    # coarse index offset + i. Coverage requires 2*i + 1 < N_f, so
    # i < (N_f - 1) / 2. Safe range: i in [0, (N_f-1)//2 - 1].
    n_covered = (N_f - 1) // 2
    for ii in range(n_covered):
        for jj in range(n_covered):
            for kk in range(n_covered):
                i = ox + ii; j = oy + jj; k = oz + kk
                if not (0 <= i < N_coarse and 0 <= j < N_coarse
                        and 0 <= k < N_coarse):
                    continue
                fi0 = 2*ii; fi1 = fi0 + 1
                fj0 = 2*jj; fj1 = fj0 + 1
                fk0 = 2*kk; fk1 = fk0 + 1
                patch = f_fine[fi0:fi1+1, fj0:fj1+1, fk0:fk1+1]
                out[i,j,k] = patch.mean()
    return out


__all__ = ["prolongate", "restrict"]


def restrict_max(f_fine, N_coarse, offset, blend=0.5):
    """Peak-preserving restriction: blend of average and max.

    When the fine patch contains a sharp peak (e.g., a compact
    Gaussian), simple averaging loses ~35% per step. This variant
    blends the arithmetic average with the maximum of the 2x2x2
    fine cells:
        out = (1 - blend) * avg + blend * max

    blend=0 -> pure averaging (matches restrict)
    blend=1 -> pure max (peak-preserving, but can bias smooth regions)
    blend=0.5 -> compromise (default)

    Returns NaN outside the fine region, same convention as restrict.
    """
    N_f = f_fine.shape[0]
    ox, oy, oz = offset
    out = np.full((N_coarse, N_coarse, N_coarse), np.nan)
    n_covered = (N_f - 1) // 2
    for ii in range(n_covered):
        for jj in range(n_covered):
            for kk in range(n_covered):
                i = ox + ii; j = oy + jj; k = oz + kk
                if not (0 <= i < N_coarse and 0 <= j < N_coarse
                        and 0 <= k < N_coarse):
                    continue
                fi0 = 2*ii; fi1 = fi0 + 1
                fj0 = 2*jj; fj1 = fj0 + 1
                fk0 = 2*kk; fk1 = fk0 + 1
                patch = f_fine[fi0:fi1+1, fj0:fj1+1, fk0:fk1+1]
                avg = patch.mean()
                mx = patch.max()
                out[i,j,k] = (1.0 - blend) * avg + blend * mx
    return out


__all__ += ["restrict_max"]
