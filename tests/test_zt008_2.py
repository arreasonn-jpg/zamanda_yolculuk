
"""ZT-008.2 tests."""
import numpy as np
from zt005.run_zt008_2 import (
    coil_resistance, joule_heating, power_breakdown, total_power,
    battery_runtime, thermal_steady_state, required_cooling_coeff,
    heat_flux_density, AMBIENT_T, COOLING_TARGET_T,
)


def test_coil_resistance_reasonable():
    R = coil_resistance()
    assert 0.01 < R < 1.0


def test_joule_heating_positive():
    P = joule_heating()
    assert P > 0


def test_power_breakdown_complete():
    b = power_breakdown()
    for key in ["coil_joule", "driver_loss", "control", "safety"]:
        assert key in b
        assert b[key] >= 0


def test_total_power_matches_sum():
    b = power_breakdown()
    assert abs(total_power() - sum(b.values())) < 1e-9


def test_battery_runtime_is_positive():
    r = battery_runtime(capacity_wh=100, power_w=50)
    assert r > 0


def test_battery_runtime_scales():
    r1 = battery_runtime(capacity_wh=100, power_w=50)
    r2 = battery_runtime(capacity_wh=100, power_w=100)
    assert abs(r1 / r2 - 2.0) < 1e-9


def test_thermal_steady_state_above_ambient():
    T = thermal_steady_state(10.0, cooling_coeff=1.0)
    assert T > AMBIENT_T


def test_required_cooling_for_target():
    hA = required_cooling_coeff(10.0, COOLING_TARGET_T)
    assert hA > 0


def test_heat_flux_positive():
    q = heat_flux_density(10.0)
    assert q > 0
