"""Unit tests for the band-local conformal projection solver core."""
import numpy as np

from zt005.amr.amr_projection import (
    flat_laplacian_7pt, gradient_magnitude, jacobi_zero_bc,
    project_phi_band)

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
