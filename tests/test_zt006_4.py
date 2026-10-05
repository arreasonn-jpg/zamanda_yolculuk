
"""ZT-006.4 tests."""
import numpy as np
from zt005.run_zt006_4 import (
    required_amplification, current_for_1s_per_year,
    scaling_table, comparison_benchmarks,
)

SECONDS_PER_YEAR = 3.15576e7


def test_amplification_factor_is_huge():
    """From 1e-48 to 6e-8, we need ~1e20 current amplification."""
    base_h00 = 1e-48
    amp = required_amplification(base_h00)
    assert 1e19 < amp < 1e21


def test_required_current_is_absurd():
    h_base = 6e-48
    h_t, I_req = current_for_1s_per_year(h_base)
    assert 1e21 < I_req < 1e23


def test_scaling_table_monotone():
    rows = scaling_table(1e-48, "test")
    deltas = [r["delta_tau_per_year_s"] for r in rows]
    assert all(deltas[i] < deltas[i+1] for i in range(len(deltas)-1))


def test_target_1s_per_year_reachable_only_with_1e20():
    """Requires I ~ 1e22 A, which is >> magnetar currents."""
    h_base = 6e-48
    _, I_req = current_for_1s_per_year(h_base)
    assert I_req > 1e21


def test_benchmarks_include_gps_and_target():
    b = comparison_benchmarks()
    assert "GPS satellites (net)" in b
    assert "Target (1 s / year)" in b
    assert b["Target (1 s / year)"] == 1.0


def test_gps_dilation_is_microseconds_per_day():
    b = comparison_benchmarks()
    gps_per_day = b["GPS satellites (net)"] / 365.25
    # 38.6 us/day
    assert 3e-5 < gps_per_day < 5e-5
