import numpy as np
import pytest
from zt005.quantum_inequality import (
    sampling_function, ford_roman_bound, quantum_inequality_value,
    check_ford_roman,
)


def test_sampling_normalized():
    """Lorentzian integrates to 1 over R."""
    t = np.linspace(-1e4, 1e4, 100001)
    tau = 1.0
    g = sampling_function(t, tau)
    assert np.trapezoid(g, t) == pytest.approx(1.0, abs=1e-3)


def test_ford_roman_bound_negative():
    assert ford_roman_bound(1.0) < 0


def test_ford_roman_bound_scaling():
    """Bound ~ tau^-4."""
    b1 = ford_roman_bound(1.0)
    b2 = ford_roman_bound(2.0)
    assert b2 / b1 == pytest.approx(1.0/16.0, rel=1e-6)


def test_positive_rho_no_violation():
    """All positive rho: no violation."""
    t = np.linspace(-10, 10, 201)
    rho = np.ones_like(t)
    result = check_ford_roman(rho, t, tau=1.0)
    assert result["violates"] is False
    assert result["value"] > 0


def test_mild_negative_no_violation():
    """Small negative dip: value > bound."""
    t = np.linspace(-10, 10, 201)
    rho = -1e-10 * np.ones_like(t)
    result = check_ford_roman(rho, t, tau=1.0)
    # bound = -3/(32 pi^2) ~ -0.0095. -1e-10 > -0.0095 -> no violation
    assert result["violates"] is False


def test_severe_negative_violates():
    """Large sustained negative rho -> violation."""
    t = np.linspace(-10, 10, 201)
    rho = -100.0 * np.ones_like(t)
    result = check_ford_roman(rho, t, tau=1.0)
    assert result["violates"] is True


def test_bound_scales_with_field_count():
    t = np.linspace(-10, 10, 201)
    rho = -0.01 * np.ones_like(t)
    r1 = check_ford_roman(rho, t, tau=1.0, N_field_components=1)
    r2 = check_ford_roman(rho, t, tau=1.0, N_field_components=2)
    # Same value, doubled bound
    assert r2["bound"] == pytest.approx(2.0 * r1["bound"])
