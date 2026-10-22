"""Band-local conformal projection solve (v13.2-alpha1).

Goal: after coarse-to-fine band injection, correct phi inside the
injected band so the fine Hamiltonian constraint is reduced, leaving
all other fields untouched.

Identity used (conformal geometric chain):
    H = e^{-4 phi} (Rtilde - 8 Dcaret^2 phi - 8 |nabla phi|^2)
        + K^2 - At_ij At^ij - 16 pi rho
so to first order the second-derivative response of H to phi is
-dH/dphi ~ 8 e^{-4 phi} Lap(phi) in the flat-Laplacian approximation
(with gamma-tilde Christoffel terms lagged; see the plan doc for the
exact covariant mapping, deferred to a later increment).

Guard: the quadratic |nabla phi|^2 term is large near punctures, so
cells whose gradient magnitude exceeds `mask_threshold` are excluded
from the source (their phi stays as prolongated).

Solve core: damped point-Jacobi on the 7-point Laplacian with
Dirichlet (zero) boundary, which is cheap on a 17^3 band.

Honesty: this is a FIRST-ORDER, masked, flat-Laplacian approximation.
It may reduce H in the band, be neutral, or make it worse — all
outcomes are recorded (see docs/AMR_PROJECTED_INJECTION_PLAN.md).
"""
from __future__ import annotations
import numpy as np


def flat_laplacian_7pt(u, h):
    """7-point central Laplacian on a cubic grid with spacing h."""
    c = u[1:-1, 1:-1, 1:-1]
    s = (u[:-2, 1:-1, 1:-1] + u[2:, 1:-1, 1:-1]
         + u[1:-1, :-2, 1:-1] + u[1:-1, 2:, 1:-1]
         + u[1:-1, 1:-1, :-2] + u[1:-1, 1:-1, 2:])
    return (s - 6.0 * c) / (h * h)


def jacobi_zero_bc(source, h, sweeps=20, omega=0.8):
    """Solve Lap(u) = source on the interior with u = 0 on the boundary.

    Damped point-Jacobi: u <- (1-omega) u + omega (sum_neighbors
    + h^2 f) / 6. Returns u of the same shape as source.
    """
    u = np.zeros_like(source)
    n = source.shape[0]
    if n < 5:
        raise ValueError("grid too small for a 7-point stencil")
    for _ in range(sweeps):
        s = (u[:-2, 1:-1, 1:-1] + u[2:, 1:-1, 1:-1]
             + u[1:-1, :-2, 1:-1] + u[1:-1, 2:, 1:-1]
             + u[1:-1, 1:-1, :-2] + u[1:-1, 1:-1, 2:])
        u[1:-1, 1:-1, 1:-1] = ((1.0 - omega) * u[1:-1, 1:-1, 1:-1]
                               + omega * (s - h * h * source[1:-1,
                                                            1:-1,
                                                            1:-1]) / 6.0)
    return u


def gradient_magnitude(phi, h):
    """Central-difference |grad phi| with one-sided edges."""
    g = np.gradient(phi, h)
    return np.sqrt(g[0] ** 2 + g[1] ** 2 + g[2] ** 2)


def project_phi_band(state_phi, h_f, hamiltonian_res, n_layers=4,
                     sweeps=20, omega=0.8, mask_threshold=0.5):
    """One projection step for the phi band.

    Source is defined only on in-band cells (n_layers deep from each
    face, additionally masked by a gradient-magnitude guard); the
    Jacobi solve runs on the full grid interior with u = 0 on the
    boundary, and only the in-band part of delta_phi is applied.

    Returns the updated phi array (copy), and the correction applied.
    """
    n_layers = min(n_layers, state_phi.shape[0] // 4)
    source = np.zeros_like(state_phi)
    safe = gradient_magnitude(state_phi, h_f) < mask_threshold
    body = safe.copy()
    for k in range(n_layers):
        body[k, ...] = False
        body[-1 - k, ...] = False
        body[:, k, ...] = False
        body[:, -1 - k, ...] = False
        body[:, :, k, ...] = False
        body[:, :, -1 - k, ...] = False
    band = safe & ~body                # in-band cells only
    src = np.where(band, np.exp(4.0 * state_phi) * hamiltonian_res / 8.0,
                   0.0)
    dphi = jacobi_zero_bc(src, h_f, sweeps=sweeps, omega=omega)
    out = state_phi.copy()
    out[band] = state_phi[band] + dphi[band]
    return out, float(np.max(np.abs(dphi[band]))) if band.any() else 0.0
