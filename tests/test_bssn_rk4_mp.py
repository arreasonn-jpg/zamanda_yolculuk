import numpy as np
from zt005.bssn_state_mp import minkowski_mp_state
from zt005.bssn_rk4_mp import rk4_step_mp


def test_minkowski_stationary():
    s = minkowski_mp_state(5)
    s2 = rk4_step_mp(s, h=0.1, dt=1e-4)
    assert np.allclose(s2.phi, 0.0, atol=1e-12)
    assert np.allclose(s2.B, 0.0, atol=1e-12)
    assert np.allclose(s2.beta, 0.0, atol=1e-12)


def test_returns_B():
    s = minkowski_mp_state(5)
    s2 = rk4_step_mp(s, h=0.1, dt=1e-4)
    assert s2.B.shape == (5,5,5,3)


def test_small_perturbation_finite():
    N = 7
    s = minkowski_mp_state(N)
    x = np.linspace(-1, 1, N)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    s.phi = 1e-4 * np.exp(-(X**2+Y**2+Z**2)/0.3)
    s2 = rk4_step_mp(s, h=2.0/(N-1), dt=1e-5)
    assert np.all(np.isfinite(s2.phi))
    assert np.all(np.isfinite(s2.B))
