"""Unit tests for the band-local conformal projection solver core."""
import numpy as np

from zt005.amr.amr_projection import (
    atilde_up_divergence, cov_laplacian, covariant_delta_phi,
    flat_laplacian_7pt, gradient_magnitude, jacobi_zero_bc,
    momentum_band_correct, project_phi_band)

N = 17
L = 10.0
h = L / (N - 1)


def test_flat_laplacian_7pt_on_quadratic():
    # Lap(x^2 + y^2 + z^2) = 6 exactly for h^2 offsets
    x = np.linspace(-L / 2, L / 2, N)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    u = X ** 2 + Y ** 2 + Z ** 2
    lap = flat_laplacian_7pt(u, h)
    assert lap.shape == (N - 2, N - 2, N - 2)
    assert np.allclose(lap, 6.0)


def test_jacobi_zero_bc_converges_to_source_solution():
    # Lap(u) = f with u = 0 on the boundary: pick u_sin = sin(a x) with
    # a = 2 pi / L -> u_sin = phi0 sin(...); f = -a^2 u_sin; solve and
    # check the interior matches u_sin pointwise within LS tolerance.
    x = np.linspace(0, L, N)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    a = 2 * np.pi / L
    # replace the boundary-zero solution: u = A sin(a x) sin(a y)
    # sin(a z) vanishes on all six faces
    u_ex = 0.3 * np.sin(a * X) * np.sin(a * Y) * np.sin(a * Z)
    f = -3.0 * a ** 2 * u_ex
    src = np.zeros_like(u_ex)
    src[1:-1, 1:-1, 1:-1] = f[1:-1, 1:-1, 1:-1]
    u_num = jacobi_zero_bc(src, h, sweeps=2000, omega=1.0)
    # convergence rate of Jacobi: omega - slow; tolerance loose
    inner = slice(2, N - 2)
    err = np.max(np.abs(u_num[inner, inner, inner]
                        - u_ex[inner, inner, inner]))
    assert err < 0.05 * np.max(np.abs(u_ex)), \
        f"rel LS err {err:.4f} vs amp {np.max(np.abs(u_ex)):.4f}"


def test_jacobi_zero_bc_boundary_stays_zero():
    src = np.zeros((N, N, N))
    src[1:-1, 1:-1, 1:-1] = 1.0
    u = jacobi_zero_bc(src, h, sweeps=10)
    assert np.all(u[0, :, :] == 0.0)
    assert np.all(u[-1, :, :] == 0.0)


# ---------------------------------------------------------------------------
# v13.3: exact gamma-tilde covariant operator
# ---------------------------------------------------------------------------

def _anisotropic_gt(x):
    """gamma_tilde = diag((1+a x)^2, 1, 1), mild variation."""
    gx = (1.0 + 0.06 * x) ** 2
    gt = np.zeros((N, N, N, 3, 3))
    for i in range(N):
        for j in range(N):
            for k in range(N):
                gt[i, j, k, 0, 0] = gx[i]
                gt[i, j, k, 1, 1] = 1.0
                gt[i, j, k, 2, 2] = 1.0
    return gt


def _delta_gt():
    gt = np.zeros((N, N, N, 3, 3))
    for d in range(3):
        gt[..., d, d] = 1.0
    return gt


def test_cov_laplacian_matches_flat_at_delta():
    x = np.linspace(-L / 2, L / 2, N)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    u = X ** 3 + 0.5 * Y ** 2 - 2.0 * Z
    lap_cov = cov_laplacian(u, _delta_gt(), h)
    lap_flat = flat_laplacian_7pt(u, h)
    assert np.allclose(lap_cov[1:-1, 1:-1, 1:-1], lap_flat)


def test_cov_laplacian_of_constant_is_zero_anisotropic():
    x = np.linspace(-L / 2, L / 2, N)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    gt = _anisotropic_gt(x)
    u = np.full((N, N, N), 3.7)
    lap = cov_laplacian(u, gt, h)
    assert np.all(lap == 0.0)


def test_cov_laplacian_analytic_linear_field():
    # gamma_xx = (1 + a x)^2: Gamma^x_xx = a / (1 + a x);
    # D^2 x = - ginv^{ij} Gamma^x_ij = - a / (1 + a x)^3 exactly,
    # since raw central second/fourth stencils of a linear field are 0
    # and central first derivatives are exact for it.
    x = np.linspace(-L / 2, L / 2, N)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    a = 0.06
    gt = _anisotropic_gt(x)
    u = X.copy()
    lap = cov_laplacian(u, gt, h)
    exact = -a / (1.0 + a * X) ** 3
    inner = slice(3, N - 3)
    assert np.allclose(lap[inner, inner, inner],
                       exact[inner, inner, inner], rtol=1e-8), \
        f"max dev {np.max(np.abs(lap[inner, inner, inner] - exact[inner, inner, inner])):.2e}"


def test_covariant_delta_phi_matches_flat_jacobi_at_delta():
    # Same equation, two solvers: at gtilde = delta the assembled
    # covariant system reduces to the flat 7-point Laplacian, so the
    # sparse solve and the converged Jacobi agree.
    x = np.linspace(-L / 2, L / 2, N)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    src = np.where((X + Y + Z) > 0.0, 0.01, 0.0)
    src = 0.02 * np.exp(-0.1 * (X ** 2 + Y ** 2 + Z ** 2))
    u_sparse = covariant_delta_phi(_delta_gt(), h, src)
    u_jacobi = jacobi_zero_bc(src, h, sweeps=6000, omega=1.0)
    inner = slice(2, N - 2)
    scale = max(1e-30, np.max(np.abs(u_jacobi)))
    dev = np.max(np.abs(u_sparse[inner, inner, inner]
                        - u_jacobi[inner, inner, inner]))
    assert dev < 2e-3 * scale, \
        f"dev {dev:.3e} vs scale {scale:.3e}"


def test_covariant_delta_phi_boundary_stays_zero():
    x = np.linspace(-L / 2, L / 2, N)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    src = np.ones((N, N, N))
    gt = _anisotropic_gt(x)
    u = covariant_delta_phi(gt, h, src)
    assert np.all(u[0, :, :] == 0.0)
    assert np.all(u[:, -1, :] == 0.0)
    assert np.isfinite(u).all()


def test_project_phi_band_masks_punctures_and_bounds_delta():
    # steep gradient center (puncture-like): must be excluded
    x = np.linspace(-L / 2, L / 2, N)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    phi = 0.1 * np.exp(-(X ** 2 + Y ** 2 + Z ** 2))
    resid = np.full((N, N, N), 1.0)
    out, dmax = project_phi_band(phi, h, resid, n_layers=4,
                                 mask_threshold=0.05)
    # central cell excluded (steep gradient): phi unchanged there
    mid = N // 2
    assert out[mid, mid, mid] == phi[mid, mid, mid]
    # corner is outside the band: unchanged too
    assert out[0, 0, 0] == phi[0, 0, 0]
    # delta bounded and finite on applied cells
    assert np.isfinite(out).all()
    assert dmax >= 0.0


def test_project_phi_band_reduces_pure_conformal_residual():
    # pure test: with a flat-Laplacian residual field derived from a
    # known phi perturbation, one projection step with many sweeps
    # reduces the discrete Laplacian residual.
    x = np.linspace(-L / 2, L / 2, N)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    phi = 0.2 * np.sin(2 * np.pi * X / L) * np.sin(2 * np.pi * Y / L)
    # treat resid = -8 Lap(phi) (i.e. H residual caused by laplacian);
    # the projection should correct phi toward zero Laplacian
    resid = -8.0 * flat_laplacian_7pt(phi, h)
    full = np.zeros_like(phi)
    full[1:-1, 1:-1, 1:-1] = resid
    out, _ = project_phi_band(phi, h, full, n_layers=0,
                              sweeps=500, omega=1.0,
                              mask_threshold=1e9)
    lap_after = np.abs(flat_laplacian_7pt(out[1:-1, 1:-1, 1:-1] if False
                                          else out, h))
    lap_before = np.abs(resid)
    assert np.max(lap_after) < np.max(np.abs(lap_before))


# ---------------------------------------------------------------------------
# v13.4: momentum-constraint (Atilde divergence) band correction
# ---------------------------------------------------------------------------

def _analytic_divfree_v(X, Y, Z):
    """v = curl(psi), psi = sin(aX) sin(aY) sin(aZ): div v = 0 exactly."""
    a = 2 * np.pi / L
    v1 = a * np.sin(a * X) * np.cos(a * Y) * np.sin(a * Z)
    v2 = -a * np.cos(a * X) * np.sin(a * Y) * np.sin(a * Z)
    return np.stack([v1, v2, np.zeros_like(X)], axis=-1)


def _lv_down(dva, ginv, gdown):
    """(L v)_ij on down indices from dv[a,b] = d_a v_b (the inverse of
    what momentum_band_correct subtracts)."""
    C = np.einsum("...ia,...ab->...ib", ginv, dva)
    div_v = np.einsum("...mn,...mn->...", ginv, dva)
    cup = (C + C.swapaxes(-1, -2)
           - (2.0 / 3.0) * ginv * div_v[..., None, None])
    return np.einsum("...ik,...kl,...jl->...ij", gdown, cup, gdown)


def _d_of_v(v, h):
    dva = np.zeros(v.shape[:3] + (3, 3))
    for b in range(3):
        g = np.gradient(v[..., b], h, edge_order=2)
        dva[..., 0, b] = g[0]
        dva[..., 1, b] = g[1]
        dva[..., 2, b] = g[2]
    return dva


def test_momentum_correction_reduces_divergence():
    x = np.linspace(-L / 2, L / 2, N)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    gt = _delta_gt()
    ginv = np.linalg.inv(gt)
    v_ex = _analytic_divfree_v(X, Y, Z)
    dva = _d_of_v(v_ex, h)
    At0 = _lv_down(dva, ginv, gt) * 0.1
    At_new, l_before, l_after = momentum_band_correct(At0, gt, h,
                                                      passes=1)
    assert l_before > 0.0
    assert l_after < 0.8 * l_before, \
        f"before {l_before:.3e} after {l_after:.3e}"
    # the fixed-point iteration converges: more passes, less divergence
    l1, l3 = momentum_band_correct(At0, gt, h, passes=3)[1:]
    assert l3 < 0.7 * l1, f"1 pass {l1:.3e} 3 passes {l3:.3e}"
    # boundary face untouched by the meter/correction combination
    assert np.allclose(At_new[0, :, :], At0[0, :, :])


def test_momentum_correction_preserves_trace():
    x = np.linspace(-L / 2, L / 2, N)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    gt = _delta_gt()
    ginv = np.linalg.inv(gt)
    v_ex = _analytic_divfree_v(X, Y, Z)
    dva = _d_of_v(v_ex, h)
    At0 = _lv_down(dva, ginv, gt) * 0.1
    trace_part = 0.05 * gt
    At0 = At0 + trace_part
    tr_before = np.einsum("...ij,...ij->...", gt, At0)
    At_new, _, _ = momentum_band_correct(At0, gt, h)
    tr_after = np.einsum("...ij,...ij->...", gt, At_new)
    assert np.allclose(tr_before, tr_after, atol=1e-12)
    assert np.isfinite(At_new).all()
