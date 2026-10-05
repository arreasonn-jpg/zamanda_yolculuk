
"""ZT-007B.5 tests - Alcubierre."""
import numpy as np
from zt005.run_zt007_b5 import (
    alcubierre_energy, alcubierre_energy_density,
    casimir_max_density, mass_equivalent,
    C, MSUN,
)


def test_energy_negative():
    E = alcubierre_energy(0.1 * C, 10.0)
    assert E < 0


def test_energy_scales_v_squared():
    E1 = abs(alcubierre_energy(0.1 * C, 10.0))
    E2 = abs(alcubierre_energy(0.2 * C, 10.0))
    assert abs(E2 / E1 - 4.0) < 1e-10


def test_energy_scales_R_squared():
    E1 = abs(alcubierre_energy(0.1 * C, 10.0))
    E2 = abs(alcubierre_energy(0.1 * C, 20.0))
    assert abs(E2 / E1 - 4.0) < 1e-10


def test_1m_bubble_requires_many_earth_masses():
    # 1 m bubble at 0.1c: |E| ~ 9e57 J = ~ 1.7e16 Earth masses
    E = abs(alcubierre_energy(0.1 * C, 1.0))
    m = mass_equivalent(E)
    m_earth = m / 5.972e24
    assert 1e15 < m_earth < 1e18


def test_100m_superluminal_requires_solar_mass():
    # 100 m bubble at c: |E| ~ 9e63 J = ~ 5e16 Solar masses
    E = abs(alcubierre_energy(1.0 * C, 100.0))
    m = mass_equivalent(E)
    m_sun = m / MSUN
    assert 1e15 < m_sun < 1e18


def test_energy_density_negative():
    rho = alcubierre_energy_density(0.1 * C, 10.0)
    assert rho < 0


def test_casimir_gap_huge():
    """Warp bubble needs ~1e35+ more negative energy than Casimir."""
    rho_warp = abs(alcubierre_energy_density(0.1 * C, 1.0))
    rho_cas  = casimir_max_density(1e-9)
    gap = rho_warp / rho_cas
    assert gap > 1e20
