
"""ZT-007B.2 tests - Tipler cylinder."""
import numpy as np
from zt005.run_zt007_b2 import (
    ctc_threshold, critical_omega, critical_rho, critical_R,
    materials_table, C2_OVER_G,
)


def test_c2_over_G_order():
    # c^2/G ~ 1.35e27 kg/m
    assert 1e27 < C2_OVER_G < 2e27


def test_ctc_threshold_satisfied_when_product_huge():
    # Huge density and rotation
    th = ctc_threshold(1e20, 1e10, 1.0)
    assert th["satisfied"] is True
    assert th["margin"] > 1.0


def test_ctc_threshold_not_satisfied_for_water():
    # Water, 1 rad/s, R=1m - hopelessly low
    th = ctc_threshold(1e3, 1.0, 1.0)
    assert th["satisfied"] is False
    assert th["margin"] < 1e-20


def test_critical_omega_for_neutron_star_density():
    # rho = 10^18 kg/m^3, R=1 m -> omega_crit = c^2/(2G rho R^3)
    wc = critical_omega(1e18, 1.0)
    expected = (2.99792458e8)**2 / (6.67430e-11) / (2 * 1e18)
    assert abs(wc - expected) / expected < 1e-6
    # Around 6.7e8 rad/s
    assert 1e8 < wc < 1e9


def test_critical_density_inverse():
    # At omega = 10^3 rad/s, R=1 m, what density?
    rc = critical_rho(1e3, 1.0)
    assert 1e23 < rc < 1e25


def test_critical_R_neutron_star():
    # Nuc density + ms spin: R_crit should be small-ish
    Rc = critical_R(1e18, 1e3)
    # Check it's the cube root of a huge number
    assert 1e1 < Rc < 1e3


def test_materials_include_neutron_star():
    m = materials_table()
    assert "Neutron star core" in m
    assert m["Neutron star core"] > 1e17


def test_increasing_R_helps_ctc():
    # Same rho, omega, larger R -> higher margin
    m1 = ctc_threshold(1e18, 1e3, 1.0)["margin"]
    m2 = ctc_threshold(1e18, 1e3, 10.0)["margin"]
    assert m2 > m1 * 100  # R^3 scaling
