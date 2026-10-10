import numpy as np
from zt005.bssn_state import minkowski_state
from zt005.bssn_v2 import rk4_step_v2


def test_v2_minkowski_stationary():
    s = minkowski_state(5)
    s2 = rk4_step_v2(s, h=0.1, dt=1e-4)
    assert np.allclose(s2.phi, 0.0, atol=1e-12)
    assert np.allclose(s2.K, 0.0, atol=1e-12)
    assert np.allclose(s2.Gamma, 0.0, atol=1e-12)
    assert np.allclose(s2.alpha, 1.0, atol=1e-12)
    assert np.allclose(s2.beta, 0.0, atol=1e-12)


def test_v2_returns_extended_state():
    s = minkowski_state(5)
    s2 = rk4_step_v2(s, h=0.1, dt=1e-4)
    assert s2.Gamma.shape == (5,5,5,3)
    assert s2.alpha.shape == (5,5,5)
    assert s2.beta.shape == (5,5,5,3)


def test_v2_small_perturbation_finite():
    s = minkowski_state(7)
    x = np.linspace(-1, 1, 7)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    s.phi = 1e-5 * np.exp(-(X**2 + Y**2 + Z**2) / 0.1)
    s2 = rk4_step_v2(s, h=0.3, dt=1e-5)
    assert np.all(np.isfinite(s2.phi))
    assert np.all(np.isfinite(s2.alpha))
    assert np.all(np.isfinite(s2.Gamma))
