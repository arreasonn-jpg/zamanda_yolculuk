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
from ..bssn_grid import christoffel_from_grid


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


# ---------------------------------------------------------------------------
# v13.3: exact gamma-tilde covariant operator
#
# The v13.2-alpha1 solve used the flat-Laplacian stand-in for the
# leading term -8 e^{-4 phi} Dtilde^2 phi of
#     H = e^{-4 phi} (Rtilde - 8 D^2 phi - 8 |grad phi|^2) + K^2 - At^2.
# Here D^2 is the covariant Laplacian of gamma-tilde,
#     D^2 u = gtmp^{ij} (d_i d_j u - Gamma^k_ij d_k u),
# with second differences taken as RAW central differences on the
# interior (identical stencil family to flat_laplacian_7pt, so the
# two operators agree exactly at gtilde = delta).
#
# Still out of scope (documented lag): the first-order-in-dphi
# cross term -16 grad(phi0) . grad(dphi) and the multiplicative
# e-fold reaction term; the guard mask keeps |grad phi| cells where
# the lag is largest out of the solve.
# ---------------------------------------------------------------------------


def _operator_coeffs(gtilde, h):
    """Batched operator coefficients for D^2 (gamma-tilde covariant).

    Returns (ginv, w): ginv[..., i, j] = (gamma_tilde^{-1})_{ij},
    w[..., k] = ginv^{ij} Gamma^k_ij (the divergence-side vector).
    """
    ginv = np.linalg.inv(gtilde)
    G, _ = christoffel_from_grid(gtilde, h)
    w = np.einsum("...ij,...kij->...k", ginv, G)
    return ginv, w


def cov_laplacian(u, gtilde, h):
    """D^2 u on the interior, zero elsewhere (raw central stencil)."""
    ginv, w = _operator_coeffs(gtilde, h)
    n = u.shape[0]
    out = np.zeros_like(u)
    c = u[1:-1, 1:-1, 1:-1]
    s = np.zeros_like(c)
    for ax in range(3):
        sl_p = [slice(1, -1)] * 3
        sl_m = [slice(1, -1)] * 3
        sl_p[ax] = slice(2, None)
        sl_m[ax] = slice(None, -2)
        g = ginv[..., ax, ax][1:-1, 1:-1, 1:-1]
        s += g * (u[tuple(sl_p)] + u[tuple(sl_m)] - 2.0 * c) / (h * h)
    pairs = [(0, 1), (0, 2), (1, 2)]
    for a, b in pairs:
        sl_pp = [slice(1, -1)] * 3
        sl_pm = [slice(1, -1)] * 3
        sl_mp = [slice(1, -1)] * 3
        sl_mm = [slice(1, -1)] * 3
        sl_pp[a] = slice(2, None)
        sl_pp[b] = slice(2, None)
        sl_pm[a] = slice(2, None)
        sl_pm[b] = slice(None, -2)
        sl_mp[a] = slice(None, -2)
        sl_mp[b] = slice(2, None)
        sl_mm[a] = slice(None, -2)
        sl_mm[b] = slice(None, -2)
        g = ginv[..., a, b][1:-1, 1:-1, 1:-1]
        s += g * (u[tuple(sl_pp)] - u[tuple(sl_pm)]
                  - u[tuple(sl_mp)] + u[tuple(sl_mm)]) / (4.0 * h * h)
    for k in range(3):
        sl_p = [slice(1, -1)] * 3
        sl_m = [slice(1, -1)] * 3
        sl_p[k] = slice(2, None)
        sl_m[k] = slice(None, -2)
        wk = w[..., k][1:-1, 1:-1, 1:-1]
        s -= wk * (u[tuple(sl_p)] - u[tuple(sl_m)]) / (2.0 * h)
    out[1:-1, 1:-1, 1:-1] = s
    return out


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


def _band_mask(phi, h, n_layers, mask_threshold):
    """In-band cells that pass the gradient guard (v13.2 semantics)."""
    n_layers = min(n_layers, phi.shape[0] // 4)
    safe = gradient_magnitude(phi, h) < mask_threshold
    body = safe.copy()
    for k in range(n_layers):
        body[k, ...] = False
        body[-1 - k, ...] = False
        body[:, k, ...] = False
        body[:, -1 - k, ...] = False
        body[:, :, k, ...] = False
        body[:, :, -1 - k, ...] = False
    return safe & ~body


def covariant_delta_phi(gtilde, h, source):
    """Solve D^2 u = source on the interior with u = 0 on the boundary.

    Variable-coefficient sparse-direct solve (scipy spsolve) over the
    full grid interior; the same raw-central stencil family as
    cov_laplacian, so at gtilde = delta the solution matches the
    flat-Laplacian solve.
    """
    import scipy.sparse as sp
    import scipy.sparse.linalg as spla

    n = source.shape[0]
    ginv, w = _operator_coeffs(gtilde, h)
    m = n - 2
    rows, cols, vals, rhs = [], [], [], []

    def lin(i, j, k):
        return (i - 1) + m * ((j - 1) + m * (k - 1))

    def add(row, i, j, k, coeff):
        if 0 < i < n - 1 and 0 < j < n - 1 and 0 < k < n - 1:
            rows.append(row)
            cols.append(lin(i, j, k))
            vals.append(coeff)

    for k in range(1, n - 1):
        for j in range(1, n - 1):
            for i in range(1, n - 1):
                r = lin(i, j, k)
                gi = ginv[i, j, k]
                wk = w[i, j, k]
                for ax in range(3):
                    off = [0, 0, 0]
                    off[ax] = 1
                    p = [i + off[0], j + off[1], k + off[2]]
                    mmm = [i - off[0], j - off[1], k - off[2]]
                    g = gi[ax, ax]
                    add(r, *p, g / (h * h))
                    add(r, *mmm, g / (h * h))
                    add(r, i, j, k, -2.0 * g / (h * h))
                for a, b in [(0, 1), (0, 2), (1, 2)]:
                    g4 = gi[a, b] / (4.0 * h * h)
                    pp = [i, j, k]
                    pm = [i, j, k]
                    mp = [i, j, k]
                    mm = [i, j, k]
                    pp[a] += 1
                    pp[b] += 1
                    pm[a] += 1
                    pm[b] -= 1
                    mp[a] -= 1
                    mp[b] += 1
                    mm[a] -= 1
                    mm[b] -= 1
                    add(r, *pp, g4)
                    add(r, *pm, -g4)
                    add(r, *mp, -g4)
                    add(r, *mm, g4)
                for kk in range(3):
                    off = [0, 0, 0]
                    off[kk] = 1
                    p = [i + off[0], j + off[1], k + off[2]]
                    mmm = [i - off[0], j - off[1], k - off[2]]
                    add(r, *p, -wk[kk] / (2.0 * h))
                    add(r, *mmm, wk[kk] / (2.0 * h))
                rhs.append(source[i, j, k])
    A = sp.csr_matrix((vals, (rows, cols)),
                      shape=(m * m * m, m * m * m))
    u_int = spla.spsolve(A, np.array(rhs))
    u = np.zeros_like(source)
    for k in range(1, n - 1):
        for j in range(1, n - 1):
            for i in range(1, n - 1):
                u[i, j, k] = u_int[lin(i, j, k)]
    return u


def project_phi_band_covariant(state_phi, gtilde, h_f, hamiltonian_res,
                               n_layers=4, mask_threshold=0.5):
    """v13.3: one projection step with the exact covariant operator.

    Same band/guard semantics as project_phi_band, but the solve is
    the assembled variable-coefficient D^2 system (gamma-tilde
    consistent instead of the flat-Laplacian stand-in). Returns the
    updated phi (copy) and max |dphi| applied on the band.
    """
    band = _band_mask(state_phi, h_f, n_layers, mask_threshold)
    src = np.where(band, np.exp(4.0 * state_phi) * hamiltonian_res / 8.0,
                   0.0)
    dphi = covariant_delta_phi(gtilde, h_f, src)
    out = state_phi.copy()
    out[band] = state_phi[band] + dphi[band]
    return out, float(np.max(np.abs(dphi[band]))) if band.any() else 0.0
