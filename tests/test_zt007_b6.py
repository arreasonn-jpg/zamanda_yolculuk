
"""ZT-007B.6 tests - Godel metric."""
import numpy as np
from zt005.run_zt007_b6 import (
    godel_density, godel_ctc_radius, godel_omega_for_ctc_at_R,
    mass_of_uniform_sphere, C, G,
)


def test_godel_density_positive():
    rho = godel_density(1e-18)
    assert rho > 0


def test_godel_ctc_radius_scales_inverse_omega():
    r1 = godel_ctc_radius(1.0)
    r2 = godel_ctc_radius(10.0)
    assert abs(r2 / r1 - 0.1) < 1e-10


def test_godel_ctc_radius_magnitude():
    """For omega=1, r_CTC ~ 0.62 c ~ 1.87e8 m."""
    r = godel_ctc_radius(1.0)
    expected = 0.6232 * C
    assert abs(r - expected) / expected < 0.01


def test_inverse_roundtrip():
    omega = 1e-6
    r = godel_ctc_radius(omega)
    omega_back = godel_omega_for_ctc_at_R(r)
    assert abs(omega_back - omega) / omega < 1e-10


def test_lab_ctc_requires_huge_omega():
    """For 1 m radius, omega ~ 1.87e8 rad/s."""
    omega = godel_omega_for_ctc_at_R(1.0)
    assert 1e7 < omega < 1e10


def test_lab_ctc_density_is_astronomical():
    """rho required for 1 m CTC is ~ 8e25 kg/m^3."""
    omega = godel_omega_for_ctc_at_R(1.0)
    rho = godel_density(omega)
    assert 1e25 < rho < 1e27


def test_lab_ctc_mass_is_huge():
    """For 1 m, mass required ~ 3.5e26 kg (58 Earths)."""
    omega = godel_omega_for_ctc_at_R(1.0)
    rho = godel_density(omega)
    m = mass_of_uniform_sphere(rho, 1.0)
    assert 1e25 < m < 1e28
