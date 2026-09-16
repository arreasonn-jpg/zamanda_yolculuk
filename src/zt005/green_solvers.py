"""Independent Dual Green Solvers for Gravitational Metric Perturbation.

Addresses Section 1.2 and 1.3:
- Solver A: Tensor-product Gauss-cylindrical quadrature with singular cell regularization.
- Solver B: Independent Cartesian cell-integrated quadrature with exact boundary-face
  divergence-theorem integration for the singular observer cell:
      int_V (1/r) dV = 0.5 * oint_{d V} (r . n / r) dA.
- Cross-validation: |h_A - h_B| / |h_ref| < epsilon.
"""

from __future__ import annotations

import numpy as np

C = 299_792_458.0
G = 6.67430e-11
EINSTEIN_GREEN_FACTOR = 4.0 * G / C**4
ETA = np.diag([-1.0, 1.0, 1.0, 1.0])


def inverse_trace_reverse(hbar):
    """Convert trace-reversed hbar_munu back to physical metric perturbation h_munu."""
    h = np.asarray(hbar, float)
    tr = np.einsum("ab,...ab->...", ETA, h)
    return h - 0.5 * ETA * tr[..., None, None] if h.ndim > 2 else h - 0.5 * ETA * tr


def box_singular_integral(dx: float, dy: float, dz: float, nq: int = 12) -> float:
    """Exact integral of 1/|r| over a rectangular box [-dx/2, dx/2] x [-dy/2, dy/2] x [-dz/2, dz/2].

    Uses Gauss's divergence theorem on the vector field v = r / (2r):
        div(v) = 1/r
        int_V (1/r) dV = 0.5 * oint_{d V} (r . n / r) dA.
    Because the origin is at the box center, r is strictly bounded away from zero
    on the six planar faces, so 2D Gauss-Legendre quadrature converges with exponential accuracy.
    """
    faces = [
        (1.0, dx / 2.0, dy, dz),   # +x face
        (1.0, dx / 2.0, dy, dz),   # -x face (symmetric)
        (1.0, dy / 2.0, dx, dz),   # +y face
        (1.0, dy / 2.0, dx, dz),   # -y face
        (1.0, dz / 2.0, dx, dy),   # +z face
        (1.0, dz / 2.0, dx, dy),   # -z face
    ]
    u, w = np.polynomial.legendre.leggauss(nq)
    total = 0.0

    # 3 pairs of opposite identical faces
    for _, dist, w1, w2 in faces[::2]:
        x1 = 0.5 * w1 * u
        w1_wt = 0.5 * w1 * w
        x2 = 0.5 * w2 * u
        w2_wt = 0.5 * w2 * w
        X1, X2 = np.meshgrid(x1, x2, indexing="ij")
        W = np.outer(w1_wt, w2_wt)
        r = np.sqrt(dist * dist + X1 * X1 + X2 * X2)
        face_int = np.sum(W * (dist / r))
        # 2 symmetric faces (+dist and -dist)
        total += 2.0 * 0.5 * face_int

    return float(total)


def solver_a_tensor_product(
    obs_points: np.ndarray,
    pts: np.ndarray,
    weights: np.ndarray,
    T: np.ndarray,
    rho_m: float | None = None,
) -> np.ndarray:
    """Solver A: Tensor-product Gauss-cylindrical quadrature with spherical regularization.

    Returns physical metric perturbation h (N_obs, 4, 4).
    """
    obs = np.atleast_2d(np.asarray(obs_points, float))
    p = np.asarray(pts, float)
    w = np.asarray(weights, float)
    vals = np.asarray(T, float)

    # Estimate default rho_m from mean cell volume if not provided
    if rho_m is None:
        mean_vol = np.mean(w)
        r_eff = (3.0 * mean_vol / (4.0 * np.pi)) ** (1.0 / 3.0)
        rho_m = float(1.2 * r_eff)

    h_out = []
    for ob in obs:
        d = np.linalg.norm(p - ob[None, :], axis=1)
        near = d < rho_m
        hbar = np.zeros((4, 4), float)

        if np.any(near):
            wn = w[near]
            tloc = np.average(vals[near], axis=0, weights=wn)
            # Exact spherical cell potential: int_{r < rho} dV/r = 2*pi*rho^2
            hbar += EINSTEIN_GREEN_FACTOR * (2.0 * np.pi * rho_m**2) * tloc

        far = ~near
        if np.any(far):
            df = np.maximum(d[far], 1e-15)
            hbar += EINSTEIN_GREEN_FACTOR * np.sum((w[far] / df)[:, None, None] * vals[far], axis=0)

        h_phys = inverse_trace_reverse(hbar)
        h_out.append(h_phys)

    return np.asarray(h_out)


def solver_b_adaptive_voxel(
    obs_points: np.ndarray,
    cyl_radius: float,
    cyl_height: float,
    source_eval_fn,
    nx: int = 24,
    ny: int = 24,
    nz: int = 24,
) -> np.ndarray:
    """Solver B: Independent Cartesian voxel quadrature with boundary fraction weighting
    and exact boundary-face singular integration.

    Parameters
    ----------
    obs_points : (N_obs, 3) coordinates
    cyl_radius : cylinder radius [m] in (x, z) plane
    cyl_height : cylinder height [m] along y axis
    source_eval_fn : callable(points (N, 3)) -> (N, 4, 4) stress-energy tensor
    nx, ny, nz : Cartesian voxel resolution across bounding box

    Returns physical metric perturbation h (N_obs, 4, 4).
    """
    obs = np.atleast_2d(np.asarray(obs_points, float))
    xs = np.linspace(-cyl_radius, cyl_radius, nx)
    ys = np.linspace(-cyl_height / 2.0, cyl_height / 2.0, ny)
    zs = np.linspace(-cyl_radius, cyl_radius, nz)

    dx = float(xs[1] - xs[0])
    dy = float(ys[1] - ys[0])
    dz = float(zs[1] - zs[0])

    XX, YY, ZZ = np.meshgrid(xs, ys, zs, indexing="ij")
    grid_pts = np.column_stack([XX.ravel(), YY.ravel(), ZZ.ravel()])

    # Exact boundary volume fraction: 3x3 sub-sampling per cell in (x, z)
    sub = np.array([-0.35, 0.0, 0.35])
    frac = np.zeros(len(grid_pts))
    for sx in sub:
        for sz in sub:
            frac += ((grid_pts[:, 0] + sx * dx)**2 + (grid_pts[:, 2] + sz * dz)**2 <= cyl_radius**2)
    frac /= 9.0

    keep = frac > 0.0
    pts_b = grid_pts[keep]
    weights_b = (frac * dx * dy * dz)[keep]

    T_b = source_eval_fn(pts_b)
    singular_pot = box_singular_integral(dx, dy, dz, nq=10)

    h_out = []
    for ob in obs:
        d = np.linalg.norm(pts_b - ob[None, :], axis=1)
        is_sing = (
            (np.abs(pts_b[:, 0] - ob[0]) <= dx / 2.0)
            & (np.abs(pts_b[:, 1] - ob[1]) <= dy / 2.0)
            & (np.abs(pts_b[:, 2] - ob[2]) <= dz / 2.0)
        )

        hbar = np.zeros((4, 4), float)
        if np.any(is_sing):
            t_cell = np.mean(T_b[is_sing], axis=0)
            hbar += EINSTEIN_GREEN_FACTOR * singular_pot * t_cell

        non_sing = ~is_sing
        if np.any(non_sing):
            d_far = np.maximum(d[non_sing], 1e-15)
            hbar += EINSTEIN_GREEN_FACTOR * np.sum(
                (weights_b[non_sing] / d_far)[:, None, None] * T_b[non_sing], axis=0
            )

        h_phys = inverse_trace_reverse(hbar)
        h_out.append(h_phys)

    return np.asarray(h_out)


def compare_green_solvers(h_a: np.ndarray, h_b: np.ndarray, eps_tol: float = 0.05) -> dict:
    """Compare outputs of Solver A and Solver B.

    Kabul: |h_A - h_B| / |h_ref| < epsilon.
    """
    ha = np.asarray(h_a, float)
    hb = np.asarray(h_b, float)

    diff = ha - hb
    ref_norm = max(float(np.linalg.norm(ha.ravel())), 1e-300)
    tensor_rel_err = float(np.linalg.norm(diff.ravel()) / ref_norm)

    h00_a = float(ha[..., 0, 0].ravel()[0]) if ha.ndim >= 2 else float(ha[0, 0])
    h00_b = float(hb[..., 0, 0].ravel()[0]) if hb.ndim >= 2 else float(hb[0, 0])
    h00_ref = max(abs(h00_a), 1e-300)
    h00_rel_err = float(abs(h00_a - h00_b) / h00_ref)
    h00_abs_err = float(abs(h00_a - h00_b))

    h0i_a = ha[..., 0, 1:]
    h0i_b = hb[..., 0, 1:]
    h0i_abs_err = float(np.linalg.norm(h0i_a - h0i_b))
    h0i_ref = max(float(np.linalg.norm(h0i_a)), 1e-300)
    h0i_rel_err = float(h0i_abs_err / h0i_ref)

    hij_a = ha[..., 1:, 1:]
    hij_b = hb[..., 1:, 1:]
    hij_abs_err = float(np.linalg.norm(hij_a - hij_b))
    hij_ref = max(float(np.linalg.norm(hij_a)), 1e-300)
    hij_rel_err = float(hij_abs_err / hij_ref)

    is_passed = bool(tensor_rel_err < eps_tol and h00_rel_err < eps_tol)

    return {
        "pass": is_passed,
        "eps_tol": eps_tol,
        "tensor_rel_err": tensor_rel_err,
        "h00_rel_err": h00_rel_err,
        "h00_abs_err": h00_abs_err,
        "h0i_abs_err": h0i_abs_err,
        "h0i_rel_err": h0i_rel_err,
        "hij_rel_err": hij_rel_err,
        "h_A_frobenius": float(np.linalg.norm(ha.ravel())),
        "h_B_frobenius": float(np.linalg.norm(hb.ravel())),
    }


def uniform_sphere_exact_potential(r: float | np.ndarray, radius_m: float, rho_val: float = 1.0) -> np.ndarray:
    """Exact Newtonian gravitational potential of a uniform solid sphere."""
    r_arr = np.asarray(r, float)
    inside = r_arr <= radius_m
    pot = np.zeros_like(r_arr)
    # Inside: 2*pi*rho*(R^2 - r^2/3)
    pot[inside] = 2.0 * np.pi * rho_val * (radius_m**2 - (r_arr[inside]**2) / 3.0)
    # Outside: 4*pi*rho*R^3 / (3*r)
    pot[~inside] = (4.0 * np.pi * rho_val * radius_m**3) / (3.0 * r_arr[~inside])
    return pot
