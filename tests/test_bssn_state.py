import numpy as np
from zt005.bssn_state import BSSNState, minkowski_state


def test_minkowski_shapes():
    s = minkowski_state(5)
    assert s.phi.shape == (5,5,5)
    assert s.gtilde.shape == (5,5,5,3,3)
    assert s.K.shape == (5,5,5)
    assert s.Atilde.shape == (5,5,5,3,3)
    assert s.Gamma.shape == (5,5,5,3)
    assert s.alpha.shape == (5,5,5)
    assert s.beta.shape == (5,5,5,3)


def test_minkowski_values():
    s = minkowski_state(5)
    assert np.all(s.phi == 0.0)
    assert np.all(s.K == 0.0)
    assert np.all(s.Atilde == 0.0)
    assert np.all(s.Gamma == 0.0)
    assert np.all(s.alpha == 1.0)
    assert np.all(s.beta == 0.0)
    assert np.allclose(s.gtilde, np.eye(3))


def test_copy_independent():
    s = minkowski_state(3)
    s2 = s.copy()
    s2.phi[0,0,0] = 1.0
    assert s.phi[0,0,0] == 0.0
    assert s2.phi[0,0,0] == 1.0


def test_shape_helper():
    s = minkowski_state(7)
    assert s.shape() == (7,7,7)
