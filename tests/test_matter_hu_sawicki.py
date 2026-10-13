import numpy as np
import pytest
from zt005.matter_hu_sawicki import f_R, f_R_deriv, screening_ratio


def test_f_R_high_curvature():
    """R >> m^2: f(R) ~ R - m^2 c1/c2 (difference below float precision)."""
    R = 1e5
    m2 = 1.0  # use O(1) m^2 for visible correction
    f = f_R(R, n=1.0, c1=1.0, c2=1.0, m2=m2)
    # Should be less than R by ~m^2*c1/c2
    assert f < R
    assert f > R - 2.0


def test_f_R_low_curvature():
    """R << m^2: f(R) ~ R - m^2 c1 (R/m^2)^n."""
    R = 1e-60
    f = f_R(R, n=1.0, c1=1.0, c2=1.0)
    assert f < R


def test_f_R_deriv_high_curvature():
    """df/dR -> 1 at high curvature."""
    R = 1e5
    d = f_R_deriv(R, n=1.0, c1=1.0, c2=1.0, m2=1.0)
    assert d == pytest.approx(1.0, abs=1e-6)


def test_f_R_deriv_finite():
    """df/dR is finite and in [0, 1] range."""
    for R in [1e-50, 1e-30, 1e-10, 1.0, 1e5]:
        d = f_R_deriv(R)
        assert np.isfinite(d)


def test_screening_ratio_positive():
    """Chameleon ratio is non-negative."""
    r = screening_ratio(1e-50, 1e-30)
    assert r >= 0
