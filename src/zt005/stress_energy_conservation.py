"""Stress-energy conservation diagnostic:  partial_mu T^{mu nu} = 0.

ZT-006.3 critical item: the Einstein source must pass a conservation check
BEFORE its metric perturbation is trusted.  This module computes the
numerical four-divergence of a sampled stress-energy field on a regular
spacetime grid and reports the residual.

Coordinates x^mu = (c t, x, y, z), signature (-, +, +, +):

    partial_mu T^{mu nu} = (1/c) d/dt T^{0 nu}
                           + d/dx T^{1 nu} + d/dy T^{2 nu} + d/dz T^{3 nu}

Second-order central differences are used in all four directions; interior
points only.
"""

from __future__ import annotations

import numpy as np

C = 299_792_458.0


def divergence_residual(T_grid, dt_s, dx_m):
    """Compute the residual partial_mu T^{mu nu} on the interior grid.

    Parameters
    ----------
    T_grid : (Nt, Nx, Ny, Nz, 4, 4) sampled stress-energy tensor
    dt_s   : time sampling interval [s]
    dx_m   : spatial sampling interval [m] (isotropic grid)

    Returns dict with the residual field and summary norms.
    """
    T = np.asarray(T_grid, float)
    if T.ndim != 6 or T.shape[4:] != (4, 4):
        raise ValueError("T_grid must have shape (Nt, Nx, Ny, Nz, 4, 4)")

    # R^nu = partial_mu T^{mu nu}: each directional derivative acts on the
    # matching tensor ROW (mu = 0 for t, 1 for x, 2 for y, 3 for z).
    res = np.zeros(T.shape[:4] + (4,))
    res[1:-1, ...] += (T[2:, ..., 0, :] - T[:-2, ..., 0, :]) \
        / (2.0 * dt_s * C)
    res[:, 1:-1, ...] += (T[:, 2:, ..., 1, :] - T[:, :-2, ..., 1, :]) \
        / (2.0 * dx_m)
    res[:, :, 1:-1, ...] += (T[:, :, 2:, ..., 2, :] - T[:, :, :-2, ..., 2, :]) \
        / (2.0 * dx_m)
    res[:, :, :, 1:-1, ...] += (T[:, :, :, 2:, ..., 3, :]
                                - T[:, :, :, :-2, ..., 3, :]) / (2.0 * dx_m)

    interior = res[1:-1, 1:-1, 1:-1, 1:-1, :]
    scale = max(float(np.max(np.abs(T))), 1e-300)
    per_nu = [float(np.max(np.abs(interior[..., nu]))) for nu in range(4)]
    return {
        "residual": interior,
        "max_abs": float(np.max(np.abs(interior))),
        "max_rel": float(np.max(np.abs(interior)) / scale),
        "per_nu_abs": per_nu,
        "scale": scale,
    }


def sample_source_on_grid(source_fn, t_grid, x_grid, y_grid, z_grid):
    """Evaluate source_fn(t, x) on the tensor product grid.

    source_fn(t_array, points (N,3)) -> (N, 4, 4)
    Returns (Nt, Nx, Ny, Nz, 4, 4).
    """
    t_grid = np.asarray(t_grid, float)
    X, Y, Z = np.meshgrid(x_grid, y_grid, z_grid, indexing="ij")
    pts = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=1)
    out = np.zeros((len(t_grid), len(x_grid), len(y_grid), len(z_grid), 4, 4))
    for i, t in enumerate(t_grid):
        out[i] = source_fn(np.full(len(pts), t), pts).reshape(
            len(x_grid), len(y_grid), len(z_grid), 4, 4)
    return out
