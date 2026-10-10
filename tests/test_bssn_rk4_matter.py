import numpy as np
from zt005.bssn_state import minkowski_state
from zt005.bssn_rk4_matter import rk4_step_matter


def test_minkowski_vacuum_stationary():
    s = minkowski_state(5)
    phi_f = np.zeros((5,5,5))
    Pi_f = np.zeros((5,5,5))
    s2, phi2, Pi2 = rk4_step_matter(s, phi_f, Pi_f, h=0.1, dt=1e-4)
    assert np.allclose(s2.phi, 0.0, atol=1e-12)
    assert np.allclose(phi2, 0.0, atol=1e-12)
    assert np.allclose(Pi2, 0.0, atol=1e-12)


def test_homogeneous_pi_evolves():
    """Homogeneous Pi: d_t phi = alpha Pi -> phi grows."""
    s = minkowski_state(5)
    phi_f = np.zeros((5,5,5))
    Pi_f = np.ones((5,5,5)) * 0.01
    s2, phi2, Pi2 = rk4_step_matter(s, phi_f, Pi_f, h=0.1, dt=1e-3)
    # phi ~ dt * Pi = 1e-5
    assert np.allclose(phi2, 1e-5, atol=1e-7)
    # Pi unchanged (K=0, no gradient)
    assert np.allclose(Pi2, 0.01, atol=1e-9)


def test_shapes_preserved():
    s = minkowski_state(5)
    phi_f = np.zeros((5,5,5))
    Pi_f = np.zeros((5,5,5))
    s2, phi2, Pi2 = rk4_step_matter(s, phi_f, Pi_f, h=0.1, dt=1e-4)
    assert phi2.shape == (5,5,5)
    assert Pi2.shape == (5,5,5)
    assert s2.phi.shape == (5,5,5)
