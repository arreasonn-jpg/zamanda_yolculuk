import numpy as np
import pytest
from zt005.binary_orbit_ic import (
    circular_orbit_velocity, binary_orbit_positions,
    gaussian_binary_ic, expected_gw_frequency,
)


def test_orbit_velocity_positive():
    v = circular_orbit_velocity(1.0, 4.0)
    assert v > 0


def test_orbit_positions_symmetric():
    pos = binary_orbit_positions(4.0, 0.0, t=0.0)
    assert pos["body1"] == pytest.approx((2.0, 0.0, 0.0))
    assert pos["body2"] == pytest.approx((-2.0, 0.0, 0.0))


def test_orbit_positions_quarter_period():
    """At t = pi/(2 omega): rotation by 90 degrees."""
    omega = 1.0
    pos = binary_orbit_positions(4.0, omega, t=np.pi/2)
    assert pos["body1"][0] == pytest.approx(0.0, abs=1e-12)
    assert pos["body1"][1] == pytest.approx(2.0, abs=1e-12)


def test_gaussian_ic_shapes():
    ic = gaussian_binary_ic(N=9, L=20.0, A_amp=1e-3, sigma=1.5,
                             separation=4.0, omega=0.1)
    assert ic["phi_sc"].shape == (9, 9, 9)
    assert ic["Pi_sc"].shape == (9, 9, 9)


def test_gaussian_ic_two_peaks():
    ic = gaussian_binary_ic(N=11, L=20.0, A_amp=1e-3, sigma=1.0,
                             separation=6.0, omega=0.1)
    # Two peaks should exist (nonzero at two distinct x)
    phi = ic["phi_sc"]
    # Sum along y, z at y=z=0 slice
    mid = 11 // 2
    x_line = phi[:, mid, mid]
    # Should have two local maxima
    peaks = np.where((x_line[1:-1] > x_line[:-2]) &
                     (x_line[1:-1] > x_line[2:]))[0] + 1
    assert len(peaks) == 2


def test_gw_frequency_scaling():
    """f_GW = omega/pi."""
    assert expected_gw_frequency(np.pi) == pytest.approx(1.0)
    assert expected_gw_frequency(2*np.pi) == pytest.approx(2.0)
