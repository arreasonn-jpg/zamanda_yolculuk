import numpy as np
import pytest
from zt005.bssn_gauge import dt_alpha_1pluslog, dt_beta_gamma_driver


def test_alpha_zero_when_K_zero():
    assert dt_alpha_1pluslog(alpha=1.0, K=0.0) == 0.0


def test_alpha_decreases_for_positive_K():
    da = dt_alpha_1pluslog(alpha=1.0, K=0.1)
    assert da == pytest.approx(-0.2)


def test_beta_zero_when_Gamma_dot_and_beta_zero():
    out = dt_beta_gamma_driver(np.zeros(3), np.zeros(3), eta=1.0)
    assert np.allclose(out, 0.0)


def test_beta_driver_responds_to_Gamma_dot():
    out = dt_beta_gamma_driver(np.array([1.0, 0.0, 0.0]), np.zeros(3), eta=1.0)
    assert out[0] == pytest.approx(0.75)


def test_beta_driver_damps_beta():
    out = dt_beta_gamma_driver(np.zeros(3), np.array([1.0, 0.0, 0.0]), eta=1.0)
    assert out[0] == pytest.approx(-0.75)
