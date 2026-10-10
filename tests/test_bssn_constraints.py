import numpy as np
from zt005.bssn_state import minkowski_state
from zt005.bssn_constraints import hamiltonian_constraint_grid


def test_H_zero_minkowski():
    s = minkowski_state(5)
    H = hamiltonian_constraint_grid(s.phi, s.gtilde, s.K, s.Atilde, h=0.1)
    assert np.allclose(H, 0.0, atol=1e-10)


def test_H_responds_to_K():
    s = minkowski_state(5)
    s.K = np.ones((5,5,5)) * 0.1
    H = hamiltonian_constraint_grid(s.phi, s.gtilde, s.K, s.Atilde, h=0.1)
    # H = K^2 = 0.01
    assert np.allclose(H, 0.01, atol=1e-9)


def test_H_shape():
    s = minkowski_state(5)
    H = hamiltonian_constraint_grid(s.phi, s.gtilde, s.K, s.Atilde, h=0.1)
    assert H.shape == (5,5,5)
