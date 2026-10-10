import numpy as np
from zt005.bssn_state import minkowski_state
from zt005.bssn_full import bssn_rhs_beta1


def test_minkowski_zero():
    s = minkowski_state(5)
    rhs = bssn_rhs_beta1(s.phi, s.gtilde, s.K, s.Atilde, s.alpha, h=0.1)
    for k in rhs:
        assert np.allclose(rhs[k], 0.0, atol=1e-10), f"{k} != 0"


def test_K_bump_dK_positive():
    """K>0 uniform: dK = alpha K^2/3 (alpha=1 const -> D^2 alpha=0)."""
    s = minkowski_state(5)
    s.K = np.ones((5,5,5)) * 0.01
    rhs = bssn_rhs_beta1(s.phi, s.gtilde, s.K, s.Atilde, s.alpha, h=0.1)
    assert np.allclose(rhs["dK"], 1e-4/3.0, atol=1e-10)


def test_shapes():
    s = minkowski_state(5)
    rhs = bssn_rhs_beta1(s.phi, s.gtilde, s.K, s.Atilde, s.alpha, h=0.1)
    assert rhs["dphi"].shape == (5,5,5)
    assert rhs["dgt"].shape == (5,5,5,3,3)
    assert rhs["dK"].shape == (5,5,5)
    assert rhs["dAt"].shape == (5,5,5,3,3)
    assert rhs["dGamma"].shape == (5,5,5,3)


def test_alpha_gradient_gives_dK_contribution():
    """Non-constant alpha -> D^2 alpha != 0 -> dK != 0."""
    N = 7
    s = minkowski_state(N)
    x = np.linspace(-1, 1, N)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    s.alpha = 1.0 + 0.1 * np.exp(-(X**2+Y**2+Z**2)/0.3)
    rhs = bssn_rhs_beta1(s.phi, s.gtilde, s.K, s.Atilde, s.alpha, h=2.0/(N-1))
    assert not np.allclose(rhs["dK"], 0.0)
    assert np.all(np.isfinite(rhs["dK"]))
