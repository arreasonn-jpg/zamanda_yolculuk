import numpy as np
from zt005.bssn_gamma import dGamma_RHS
from zt005.bssn_state import minkowski_state


def _gaussian(N, amp, sigma=0.3):
    x = np.linspace(-1, 1, N)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    return amp * np.exp(-(X**2 + Y**2 + Z**2) / (2 * sigma**2))


def _smooth_atilde(N, amplitude):
    Axx = _gaussian(N, amplitude)
    At = np.zeros((N, N, N, 3, 3))
    At[..., 0, 0] = Axx
    At[..., 1, 1] = -0.5 * Axx
    At[..., 2, 2] = -0.5 * Axx
    return At


def test_dGamma_zero_on_minkowski():
    s = minkowski_state(5)
    dG = dGamma_RHS(s.Gamma, s.gtilde, s.Atilde, s.K, s.phi, s.alpha, h=0.1)
    assert np.allclose(dG, 0.0, atol=1e-12)


def test_dGamma_finite_with_alpha_gradient():
    N = 7
    s = minkowski_state(N)
    s.alpha = 1.0 + _gaussian(N, 1e-3)
    s.Atilde = _smooth_atilde(N, 1e-3)
    dG = dGamma_RHS(s.Gamma, s.gtilde, s.Atilde, s.K, s.phi, s.alpha, h=0.1)
    assert np.all(np.isfinite(dG))
    assert dG.shape == (N, N, N, 3)


def test_dGamma_scales_with_atilde():
    """Scale Atilde (not alpha) -> dG scales linearly."""
    N = 9
    s1 = minkowski_state(N)
    s1.alpha = 1.0 + _gaussian(N, 1e-3)
    s1.Atilde = _smooth_atilde(N, 1e-3)
    s2 = minkowski_state(N)
    s2.alpha = 1.0 + _gaussian(N, 1e-3)
    s2.Atilde = _smooth_atilde(N, 2e-3)
    dG1 = dGamma_RHS(s1.Gamma, s1.gtilde, s1.Atilde, s1.K, s1.phi, s1.alpha, h=0.1)
    dG2 = dGamma_RHS(s2.Gamma, s2.gtilde, s2.Atilde, s2.K, s2.phi, s2.alpha, h=0.1)
    n1 = np.max(np.abs(dG1))
    n2 = np.max(np.abs(dG2))
    assert n1 > 1e-12
    assert abs(n2/n1 - 2.0) < 0.01
