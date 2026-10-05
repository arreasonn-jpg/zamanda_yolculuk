
"""ZT-006.5 tests - geodesic integration."""
import numpy as np
from zt005.run_zt006_5 import (
    build_metric_grid, interpolate_metric,
    christoffel_numerical, integrate_geodesic, geodesic_norm,
)


def test_metric_interpolation_returns_4x4():
    g, rv, zv = build_metric_grid("G1", Nr=20, Nz=40, n_phi=4)
    g_interp = interpolate_metric(g, rv, zv, 0.5, 0.0)
    assert g_interp.shape == (4, 4)


def test_christoffel_shape():
    g, rv, zv = build_metric_grid("G1", Nr=20, Nz=40, n_phi=4)
    Gamma = christoffel_numerical(g, rv, zv, (0.3, 0.0))
    assert Gamma.shape == (4, 4, 4)


def test_rest_geodesic_norm_is_minus_one():
    """A particle at rest should have g_uu U^u U^v ~ -1."""
    g, rv, zv = build_metric_grid("G1", Nr=20, Nz=40, n_phi=4)
    sol = integrate_geodesic(g, rv, zv,
                             x0=(0, 0, 0, 0),
                             U0=(1.0, 0, 0, 0),
                             tau_max=1e-4)
    X = sol.y[:4]; U = sol.y[4:]
    norms = geodesic_norm(g, rv, zv, X, U)
    assert np.abs(norms.mean() + 1.0) < 1e-3


def test_rest_geodesic_stays_at_origin():
    """Without velocity, spatial coordinates should stay at 0."""
    g, rv, zv = build_metric_grid("G1", Nr=20, Nz=40, n_phi=4)
    sol = integrate_geodesic(g, rv, zv,
                             x0=(0, 0, 0, 0),
                             U0=(1.0, 0, 0, 0),
                             tau_max=1e-4)
    X = sol.y[:4]
    # Spatial coords should not have moved appreciably
    assert np.max(np.abs(X[1:])) < 1e-3


def test_proper_time_close_to_coordinate_time():
    """For weak field, dt/dtau ~ 1."""
    g, rv, zv = build_metric_grid("G1", Nr=20, Nz=40, n_phi=4)
    sol = integrate_geodesic(g, rv, zv,
                             x0=(0, 0, 0, 0),
                             U0=(1.0, 0, 0, 0),
                             tau_max=1e-4)
    U_t = sol.y[4]
    assert np.allclose(U_t, 1.0, atol=1e-6)
