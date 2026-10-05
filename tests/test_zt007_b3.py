
"""ZT-007B.3 tests."""
import numpy as np
from zt005.run_zt007_b3 import (
    exotic_mass_required, casimir_energy_density,
    casimir_plate_distance_for_rho, max_casimir_density,
    wormhole_casimir_gap, G, C, HBAR,
)


def test_exotic_mass_coefficient():
    # M_exotic = c^2/G * b0 = 1.35e27 * b0
    b0 = 1.0
    M = exotic_mass_required(b0)
    expected = C**2 / G
    assert abs(M - expected) / expected < 1e-6


def test_exotic_mass_is_huge():
    # For 1 mm: M ~ 1.35e24 kg (Earth mass scale)
    M = exotic_mass_required(1e-3)
    assert 1e23 < M < 1e25


def test_casimir_density_negative():
    rho = casimir_energy_density(1e-9)
    assert rho < 0


def test_casimir_density_increases_at_smaller_d():
    rho1 = abs(casimir_energy_density(1e-8))
    rho2 = abs(casimir_energy_density(1e-9))
    # rho ~ d^-4
    assert abs(rho2 / rho1 - 1e4) / 1e4 < 1e-6


def test_casimir_inverse_roundtrip():
    d = 1e-9
    rho = casimir_energy_density(d)
    d_back = casimir_plate_distance_for_rho(rho)
    assert abs(d_back - d) / d < 1e-6


def test_max_casimir_density_nanometer():
    # At d=1nm, rho = -pi^2 hbar c / (720 d^4) ~ -4.33e8 J/m^3
    rho = max_casimir_density(1e-9)
    assert 1e8 < rho < 1e9


def test_wormhole_casimir_gap_huge():
    # 1 mm wormhole gap >> 1e30
    gap = wormhole_casimir_gap(1e-3, 1e-9)["gap"]
    assert gap > 1e30


def test_casimir_does_not_scale_with_throat():
    # SMALLER wormhole needs HIGHER density (rho ~ 1/b0^2)
    # So gap_planck > gap_macro
    gap_macro = wormhole_casimir_gap(1e-3, 1e-9)["gap"]
    gap_planck = wormhole_casimir_gap(1e-35, 1e-9)["gap"]
    assert gap_planck > gap_macro * 1e50
