import numpy as np
from zt005.bssn_state import minkowski_state
from zt005.bssn_rk4_beta2 import rk4_step_beta2


def test_minkowski_stationary():
    s = minkowski_state(5)
    s2 = rk4_step_beta2(s, h=0.1, dt=1e-4)
    assert np.allclose(s2.phi, 0.0, atol=1e-12)
    assert np.allclose(s2.beta, 0.0, atol=1e-12)
    assert np.allclose(s2.alpha, 1.0, atol=1e-12)


def test_shapes_preserved():
    s = minkowski_state(5)
    s2 = rk4_step_beta2(s, h=0.1, dt=1e-4)
    assert s2.beta.shape == (5,5,5,3)
    assert s2.Gamma.shape == (5,5,5,3)


def test_small_perturbation_finite():
    N = 7
    s = minkowski_state(N)
    x = np.linspace(-1, 1, N)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    s.phi = 1e-4 * np.exp(-(X**2+Y**2+Z**2)/0.3)
    s2 = rk4_step_beta2(s, h=2.0/(N-1), dt=1e-5)
    assert np.all(np.isfinite(s2.phi))
    assert np.all(np.isfinite(s2.beta))
    assert np.all(np.isfinite(s2.Gamma))
