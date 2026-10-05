
"""ZT-006.5.1 tests - high-current sanity."""
import numpy as np
from zt005.run_zt006_5_1 import (
    scaled_h00, build_metric_grid_scaled,
)
from zt005.run_zt006_5 import integrate_geodesic, geodesic_norm


def test_scaled_h00_quadratic():
    h_base = 1e-48
    h_2x = scaled_h00(h_base, 200.0)
    assert abs(h_2x / h_base - 4.0) < 1e-10


def test_scaled_h00_at_1e22():
    h_base = 5.9e-48
    h_new = scaled_h00(h_base, 1e22)
    # h_new = 5.9e-48 * (1e20)^2 = 5.9e-8
    assert 1e-9 < h_new < 1e-7


def test_metric_scaling_builds():
    g, rv, zv, h00 = build_metric_grid_scaled("G1", 1e22,
                                              Nr=20, Nz=40, n_phi=4)
    assert g.shape == (20, 40, 4, 4)
    assert h00.shape == (20, 40)
    assert np.max(np.abs(h00)) > 1e-9


def test_geodesic_visible_at_high_current():
    """At I=1e22, geodesic deviation should be visible."""
    g, rv, zv, _ = build_metric_grid_scaled("G1", 1e22,
                                            Nr=30, Nz=60, n_phi=4)
    sol = integrate_geodesic(
        g, rv, zv, x0=(0, 0, 0, 0), U0=(1.0, 0, 0, 0),
        tau_max=1.0, n_steps=50,
    )
    X = sol.y[:4]
    # Some visible deviation from origin
    max_spatial = np.max(np.abs(X[1:]))
    assert max_spatial > 1e-12


def test_low_current_stays_minkowski():
    """At I=100 A, deviation is below float64 epsilon."""
    g, rv, zv, h00 = build_metric_grid_scaled("G1", 100.0,
                                              Nr=20, Nz=40, n_phi=4)
    # g_tt should be exactly -1.0
    assert g[0, 20, 0, 0] == -1.0
