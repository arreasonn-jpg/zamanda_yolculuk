import numpy as np
import pytest
from zt005.bssn_gauge_mp import dt_alpha_mp, dt_beta_mp, dt_B_mp


def test_alpha_zero_when_K_zero():
    assert dt_alpha_mp(1.0, 0.0) == 0.0


def test_alpha_decreases_positive_K():
    assert dt_alpha_mp(1.0, 0.1) == pytest.approx(-0.2)


def test_beta_zero_when_B_zero():
    assert np.allclose(dt_beta_mp(np.zeros(3)), 0.0)


def test_beta_responds_to_B():
    out = dt_beta_mp(np.array([1.0, 0.0, 0.0]))
    assert out[0] == pytest.approx(0.75)


def test_B_damped_when_Gamma_zero():
    out = dt_B_mp(np.zeros(3), np.array([1.0, 0.0, 0.0]), eta=2.0)
    assert out[0] == pytest.approx(-2.0)


def test_B_grows_with_dGamma():
    out = dt_B_mp(np.array([1.0, 0.0, 0.0]), np.zeros(3))
    assert out[0] == pytest.approx(1.0)
