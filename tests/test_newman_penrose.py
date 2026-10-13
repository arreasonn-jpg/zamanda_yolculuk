import numpy as np
import pytest
from zt005.newman_penrose import (
    tt_polarizations, psi4_from_hddot, psi4_amplitude,
    gw_energy_flux, strain_from_psi4,
)


def test_tt_plus_polarization():
    """Pure + polarization: h = diag(h+, -h+, 0) in (e1, e2, k) frame."""
    h = np.diag([1.0, -1.0, 0.0])
    k = np.array([0.0, 0.0, 1.0])
    e1 = np.array([1.0, 0.0, 0.0])
    e2 = np.array([0.0, 1.0, 0.0])
    hp, hc = tt_polarizations(h, k, e1, e2)
    assert hp == pytest.approx(2.0)  # e1 h e1 - e2 h e2 = 1 - (-1) = 2
    assert hc == pytest.approx(0.0)


def test_tt_cross_polarization():
    """Pure x polarization: h_x nonzero."""
    h = np.array([[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
    k = np.array([0.0, 0.0, 1.0])
    e1 = np.array([1.0, 0.0, 0.0])
    e2 = np.array([0.0, 1.0, 0.0])
    hp, hc = tt_polarizations(h, k, e1, e2)
    assert hp == pytest.approx(0.0)
    assert hc == pytest.approx(2.0)


def test_psi4_zero_for_static():
    """Static metric -> psi_4 = 0."""
    p = psi4_from_hddot(0.0, 0.0)
    assert abs(p) < 1e-15


def test_psi4_amplitude_omega_squared():
    """|psi_4| = 1/2 omega^2 |h|."""
    assert psi4_amplitude(1.0, 2.0) == pytest.approx(2.0)
    assert psi4_amplitude(1.0, 4.0) == pytest.approx(8.0)


def test_gw_energy_flux_positive():
    F = gw_energy_flux(1e-20, 0.0)
    assert F > 0
    assert np.isfinite(F)


def test_strain_inverse():
    """strain_from_psi4(psi4_amplitude(h, omega), omega) ~ h."""
    h = 1e-20
    omega = 1e3
    psi4 = psi4_amplitude(h, omega)
    h_recovered = strain_from_psi4(psi4, omega)
    assert h_recovered == pytest.approx(h, rel=1e-12)
