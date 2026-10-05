"""Target-time kinematic mapping sanity tests (PHASE 9 placeholder)."""
import numpy as np

from zt005.target_time_mapping import (AU_M, SIDEREAL_DAY_S, YEAR_S,
                                       device_state_heliocentric,
                                       displacement_between,
                                       earth_center_heliocentric,
                                       surface_point_earth_centered)


def test_orbit_radius_is_one_au():
    for t in (0.0, 1e6, 1e7):
        r = np.linalg.norm(earth_center_heliocentric(t))
        assert abs(r - AU_M) / AU_M < 1e-12


def test_surface_point_returns_after_one_sidereal_day():
    lat, lon = 39.0, 35.0
    p0 = surface_point_earth_centered(0.0, lat, lon)
    p1 = surface_point_earth_centered(SIDEREAL_DAY_S, lat, lon)
    assert np.allclose(p0, p1, atol=1.0)   # < 1 m closure error


def test_one_year_orbit_closure():
    r0 = earth_center_heliocentric(0.0)
    r1 = earth_center_heliocentric(YEAR_S)
    assert np.allclose(r0, r1, atol=1.0)


def test_displacement_scale_over_one_day():
    # the Earth moves ~2.6e6 km along its orbit per day: any target-time
    # mapping must account for this, not just for a coordinate label
    d = displacement_between(0.0, 86400.0, lat_deg=0.0, lon_deg=0.0)
    assert 2.0e9 < d < 3.2e9
