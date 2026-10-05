
"""ZT-006.3.1 tests - precision-safe proper time."""
import numpy as np
from zt005.run_zt006_3 import (
    build_metric_from_h, proper_time_delta, compute_geodesic_diagnostics,
)


def test_metric_minkowski_when_h_zero():
    g = build_metric_from_h(np.zeros((10,10)), np.zeros((10,10)))
    assert np.allclose(g[0,0,0,0], -1.0)
    assert np.allclose(g[0,0,1,1],  1.0)


def test_proper_time_delta_zero_for_h_zero():
    d = proper_time_delta(np.zeros((5,5)))
    assert np.all(d == 0.0)


def test_proper_time_delta_recovers_small_h():
    h = np.full((5,5), 1e-48)
    d = proper_time_delta(h)
    # Should be ~ h/2 = 5e-49, NOT zero
    assert np.all(d > 0.0)
    assert np.all(np.abs(d - 5e-49) / 5e-49 < 1e-6)


def test_proper_time_delta_large_h_matches_direct():
    h = np.full((5,5), 0.5)
    d = proper_time_delta(h)
    direct = 1.0 - np.sqrt(1.0 - 0.5)
    assert np.allclose(d, direct, rtol=1e-12)


def test_diagnostics_finite_and_tiny():
    res = compute_geodesic_diagnostics("G1", 30, 60, n_phi=4)
    assert np.isfinite(res["h_00_center"])
    assert res["delta_tau_over_t_center"] > 0.0
    assert res["delta_tau_over_t_center"] < 1e-40


def test_proper_time_delta_center_equals_h00_over_two():
    """For weak field, 1-dtau/dt ~ h_00/2."""
    res = compute_geodesic_diagnostics("G1", 30, 60, n_phi=4)
    h00  = res["h_00_center"]
    frac = res["delta_tau_over_t_center"]
    assert abs(frac - h00/2.0) / (h00/2.0) < 1e-8
