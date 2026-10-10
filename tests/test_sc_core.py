import numpy as np
import pytest
from zt005.sc_core import (
    schwarzschild_radius, de_sitter_length, f_SC, rho_SC,
    ricci_scalar_SC, kretschmann_SC, junction_surface_density,
    core_parameters,
)


def test_schwarzschild_radius():
    M_sun = 1.989e30
    Rs = schwarzschild_radius(M_sun)
    assert 2.9e3 < Rs < 3.0e3  # ~2.95 km


def test_de_sitter_length_positive():
    ell = de_sitter_length(1e19)  # ~ nuclear density
    assert ell > 0
    assert np.isfinite(ell)


def test_f_SC_limits():
    M = 1.989e30
    Rs = schwarzschild_radius(M)
    ell = 1e4
    R0 = ell
    # Deep inside: f_in ~ 1 at r=0
    f_center = f_SC(np.array([1e-3]), M, R0, ell)[0]
    assert f_center == pytest.approx(1.0, abs=1e-6)
    # Far outside: f_out -> 1
    f_far = f_SC(np.array([1e10]), M, R0, ell)[0]
    assert f_far == pytest.approx(1.0, abs=1e-6)


def test_ricci_finite_at_origin():
    M = 1.989e30
    R0 = 1e4
    ell = 1e4
    r = np.array([1e-2, 1.0, 10.0, 100.0])
    R = ricci_scalar_SC(r, M, R0, ell)
    assert np.all(np.isfinite(R))


def test_kretschmann_finite_at_origin():
    M = 1.989e30
    R0 = 1e4
    ell = 1e4
    r = np.array([1e-2, 1.0, 10.0, 100.0])
    K = kretschmann_SC(r, M, R0, ell)
    assert np.all(np.isfinite(K))
    # Interior: K ~ const (de Sitter)
    assert K[0] < 1e10


def test_junction_returns_finite_for_valid_params():
    M = 1.989e30
    R0 = 1e4
    ell = 1e4
    s = junction_surface_density(M, R0, ell)
    assert np.isfinite(s)


def test_core_parameters():
    M = 1.989e30
    rho0 = 1e19
    p = core_parameters(M, rho0)
    assert p["R0"] > 0
    assert p["ell"] > 0
    assert p["Rs"] > 0
    assert "delta" in p
