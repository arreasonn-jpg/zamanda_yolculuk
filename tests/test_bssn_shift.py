import numpy as np
from zt005.bssn_state import minkowski_state
from zt005.bssn_shift import shift_terms


def test_zero_shift_zero_terms():
    s = minkowski_state(5)
    L = shift_terms(s.phi, s.gtilde, s.K, s.Atilde, s.Gamma, s.beta, h=0.1)
    for k in L:
        assert np.allclose(L[k], 0.0, atol=1e-12), f"{k} != 0"


def test_constant_shift_phi_advection():
    N = 5
    s = minkowski_state(N)
    x = np.linspace(-1, 1, N)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    s.phi = X.copy()
    s.beta = np.zeros((N,N,N,3))
    s.beta[..., 0] = 1.0
    L = shift_terms(s.phi, s.gtilde, s.K, s.Atilde, s.Gamma, s.beta, h=2.0/(N-1))
    interior = L["Lphi"][2:-2, 2:-2, 2:-2]
    assert np.allclose(interior, 1.0, atol=1e-10)


def test_shapes():
    s = minkowski_state(5)
    L = shift_terms(s.phi, s.gtilde, s.K, s.Atilde, s.Gamma, s.beta, h=0.1)
    assert L["Lphi"].shape == (5,5,5)
    assert L["Lgt"].shape == (5,5,5,3,3)
    assert L["LGamma"].shape == (5,5,5,3)
