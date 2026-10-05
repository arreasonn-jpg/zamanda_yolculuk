
"""ZT-009.3 tests."""
import numpy as np
from zt005.run_zt009_3 import (
    tier_scenarios, compute_h00_with_engineering,
    max_current_thermal, ENVELOPE, THERMAL,
    H00_TARGET, MAX_STAGES, MAX_Q,
)


def test_three_tiers():
    t = tier_scenarios()
    assert len(t) == 3


def test_tiers_increasing_power():
    t = tier_scenarios()
    currents = [x["I_A"] for x in t]
    assert currents[0] < currents[1] < currents[2]


def test_best_case_uses_superconducting():
    t = tier_scenarios()
    best = t[-1]
    assert best["I_A"] >= 1000
    assert best["eta"] > 0.95
    assert best["N_stages"] >= 5


def test_compute_h00_scales_quadratically_in_current():
    h1 = compute_h00_with_engineering(100, 0.1, 0.9)
    h2 = compute_h00_with_engineering(200, 0.1, 0.9)
    assert abs(h2 / h1 - 4.0) < 1e-9


def test_compute_h00_scales_quadratically_in_stages():
    h1 = compute_h00_with_engineering(100, 0.1, 0.9, N_stages=1)
    h2 = compute_h00_with_engineering(100, 0.1, 0.9, N_stages=10)
    assert abs(h2 / h1 - 100.0) < 1e-9


def test_compute_h00_scales_linearly_in_Q():
    h1 = compute_h00_with_engineering(100, 0.1, 0.9, Q_factor=1)
    h2 = compute_h00_with_engineering(100, 0.1, 0.9, Q_factor=1000)
    assert abs(h2 / h1 - 1000.0) < 1e-9


def test_max_current_scales_with_density():
    I1 = max_current_thermal(1.0, cooling="passive")
    I2 = max_current_thermal(1.0, cooling="cryogenic")
    assert I2 > I1 * 10


def test_best_case_still_far_from_target():
    t = tier_scenarios()
    best = t[-1]
    h = compute_h00_with_engineering(
        best["I_A"], best["R_m"], best["eta"],
        best["N_stages"], best["Q_factor"],
    )
    # Even best case must be far below target
    assert h < H00_TARGET / 1e20
