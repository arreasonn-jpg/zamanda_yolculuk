"""Lorenz-gauge consistency diagnostic:  partial_mu hbar^{mu nu} = 0.

ZT-006.3 critical item: the linearized Einstein solver is only meaningful
in a fixed gauge.  The retarded integral used in retarded_gr.py solves the
Lorenz-gauge wave equation, and the exact retarded solution of a CONSERVED
source satisfies partial_mu hbar^{mu nu} = 0 identically.  This gives a
sharp self-consistency test:

    conserved source  ->  small Lorenz residual  (solver OK)
    non-conserved source / gauge-violating ansatz -> large residual

Coordinates x^mu = (c t, x, y, z), signature (-, +, +, +), SI units.
Second-order central differences, interior points only.
"""

from __future__ import annotations

import numpy as np

C = 299_792_458.0


def lorenz_residual(hbar_grid, dt_s, dx_m):
    """Compute partial_mu hbar^{mu nu} on the interior of a spacetime grid.

    Parameters
    ----------
    hbar_grid : (Nt, Nx, Ny, Nz, 4, 4) trace-reversed perturbation

    Returns dict with the residual field and summary norms.
    """
    h = np.asarray(hbar_grid, float)
    if h.ndim != 6 or h.shape[4:] != (4, 4):
        raise ValueError("hbar_grid must have shape (Nt, Nx, Ny, Nz, 4, 4)")

    # R^nu = partial_mu hbar^{mu nu}: the directional derivative along each
    # coordinate acts on the matching tensor ROW (mu index).
    res = np.zeros(h.shape[:4] + (4,))
    res[1:-1, ...] += (h[2:, ..., 0, :] - h[:-2, ..., 0, :]) \
        / (2.0 * dt_s * C)
    res[:, 1:-1, ...] += (h[:, 2:, ..., 1, :] - h[:, :-2, ..., 1, :]) \
        / (2.0 * dx_m)
    res[:, :, 1:-1, ...] += (h[:, :, 2:, ..., 2, :] - h[:, :, :-2, ..., 2, :]) \
        / (2.0 * dx_m)
    res[:, :, :, 1:-1, ...] += (h[:, :, :, 2:, ..., 3, :]
                                - h[:, :, :, :-2, ..., 3, :]) / (2.0 * dx_m)

    interior = res[1:-1, 1:-1, 1:-1, 1:-1, :]
    scale = max(float(np.max(np.abs(h))), 1e-300)
    return {
        "residual": interior,
        "max_abs": float(np.max(np.abs(interior))),
        "max_rel": float(np.max(np.abs(interior)) / scale),
        "per_nu_abs": [float(np.max(np.abs(interior[..., nu])))
                       for nu in range(4)],
        "scale": scale,
    }
