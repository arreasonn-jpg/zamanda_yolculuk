"""Integration test for traveler geodesic evolution and causality tracking."""

import numpy as np
import pytest
from zt005.traveler_geodesic_worldline import evolve_traveler_geodesic
from zt005.benchmark_metrics import MinkowskiMetric, SchwarzschildMetric, GodelMetric
from zt005.ctc_past_loop_engine import check_one_second_backward_travel_criterion
from zt005.retarded_gr import C


def test_traveler_geodesic_flat_spacetime():
    m = MinkowskiMetric()
    x0 = np.array([0.0, 0.0, 0.0, 0.0])
    u0 = np.array([C, 0.0, 0.0, 0.0])

    res = evolve_traveler_geodesic(m, x0, u0, tau_span=2.0, n_steps=50)
    assert res["all_timelike"] is True
    assert pytest.approx(res["elapsed_coordinate_time_delta_t_s"], rel=1e-8) == 2.0
    assert pytest.approx(res["elapsed_proper_time_delta_tau_s"], rel=1e-8) == 2.0
    assert res["causality_inversion_detected"] is False
    assert res["norm_drift_rel"] < 1e-12


def test_one_second_backward_criterion():
    m = GodelMetric(omega=1.0)
    # A path reaching negative coordinate time
    curve = np.array([
        [0.0, 1.5, 0.0, 0.0],
        [C * (-0.5), 1.5, np.pi / 2, 0.0],
        [C * (-1.2), 1.5, np.pi, 0.0],
    ])
    crit = check_one_second_backward_travel_criterion(m, curve, target_dt_seconds=-1.0, periodic_coord=2)
    assert crit["reaches_past_time"] is True
    assert crit["achieved_delta_t_s"] <= -1.0
    assert crit["proper_time_positive"] is True
