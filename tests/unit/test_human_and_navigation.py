"""Unit tests for Earth navigation, Human extended model, and Energy-to-CTC scan."""

import numpy as np
import pytest
from zt005.earth_navigation_frame import (
    geodetic_to_ecef,
    ecef_to_geodetic,
    compute_arrival_mapping,
)
from zt005.human_extended_model import (
    HumanProxy,
    construct_orthonormal_tetrad,
    evaluate_human_safety_gate,
)
from zt005.energy_ctc_scan import run_energy_to_ctc_scan


def test_geodetic_ecef_roundtrip():
    lat, lon, alt = 41.0082, 28.9784, 50.0
    xyz = geodetic_to_ecef(lat, lon, alt)
    lat2, lon2, alt2 = ecef_to_geodetic(xyz)
    assert pytest.approx(lat, abs=1e-5) == lat2
    assert pytest.approx(lon, abs=1e-5) == lon2
    assert pytest.approx(alt, abs=1e-2) == alt2


def test_arrival_mapping_one_second_past():
    res = compute_arrival_mapping(initial_lat_deg=41.0, initial_lon_deg=28.0, target_dt_seconds=-1.0)
    # Earth moves ~ 29.8 km along its orbit in 1 second
    assert 29000.0 < res["earth_orbital_displacement_m"] < 31000.0
    # Earth surface rotates ~ 350 m at lat 41 in 1 second
    assert 300.0 < res["earth_spin_displacement_m"] < 400.0


def test_human_proxy_safety_gate():
    human = HumanProxy()
    assert human.mass_kg == 75.0
    assert human.height_m == 1.75

    # Safe test: zero curvature
    riemann_zero = np.zeros((4, 4, 4, 4))
    g = np.diag([-1.0, 1.0, 1.0, 1.0])
    u = np.array([1.0, 0.0, 0.0, 0.0])
    res_safe = evaluate_human_safety_gate(riemann_zero, g, u, human=human)
    assert res_safe["safety_gate_pass"] is True
    assert res_safe["head_to_toe_tidal_acceleration_g"] == 0.0

    # Lethal test: huge curvature
    riemann_lethal = np.zeros((4, 4, 4, 4))
    riemann_lethal[3, 0, 3, 0] = 1e-10  # c^2 * 1e-10 ~ 1e7 m/s^2 >> 10g
    res_lethal = evaluate_human_safety_gate(riemann_lethal, g, u, human=human)
    assert res_lethal["safety_gate_pass"] is False
    assert "FAIL" in res_lethal["verdict"]


def test_energy_ctc_scan():
    scan = run_energy_to_ctc_scan(baseline_current_a=100.0, baseline_h00=1.98e-46)
    assert scan["threshold"]["I_CTC_amperes"] > 1e24
    assert scan["threshold"]["orders_of_magnitude_energy_gap"] > 40
    # Lab current should have no CTC
    assert scan["scan_records"][0]["ctc_formed"] is False
