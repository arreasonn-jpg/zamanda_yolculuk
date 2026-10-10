import numpy as np
import pytest
from zt005.bssn_ic_kerr_full import kerr_schild_metric, kerr_schild_state


def test_null_vector_normalized():
    """k_i k^i = 1 outside ring singularity (rho > a)."""
    N = 9; L = 20.0; M = 1.0; a = 0.5
    R, H, k, alpha, bi, bu, g = kerr_schild_metric(N=N, L=L, M=M, a=a)
    x = np.linspace(-L/2, L/2, N)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    rho = np.sqrt(X**2 + Y**2)
    mask = rho > 2.0 * a  # outside ring
    knorm = np.einsum('...i,...i->...', k, k)
    assert np.allclose(knorm[mask], 1.0, atol=1e-10)


def test_H_positive():
    R, H, k, alpha, bi, bu, g = kerr_schild_metric(N=9, L=20.0, M=1.0, a=0.9)
    assert np.all(H >= 0.0)
    assert np.all(H[2:-2, 2:-2, 2:-2] > 0.0)


def test_alpha_range():
    R, H, k, alpha, bi, bu, g = kerr_schild_metric(N=9, L=20.0, M=1.0, a=0.5)
    assert np.all(alpha > 0.0)
    assert np.all(alpha <= 1.0)
    assert alpha[0,0,0] > 0.9


def test_a0_r_equals_rho():
    R, H, k, alpha, bi, bu, g = kerr_schild_metric(N=9, L=20.0, M=1.0, a=0.0)
    x = np.linspace(-10, 10, 9)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    rho = np.sqrt(X**2 + Y**2 + Z**2)
    assert np.allclose(R, rho, atol=1e-6)


def test_state_returns_valid_bssn():
    s = kerr_schild_state(N=9, L=20.0, M=1.0, a=0.5)
    assert s.phi.shape == (9,9,9)
    assert s.beta.shape == (9,9,9,3)
    assert np.all(np.isfinite(s.phi))
    assert np.all(np.isfinite(s.alpha))
    assert np.all(np.isfinite(s.K))


def test_a_too_large_raises():
    with pytest.raises(ValueError):
        kerr_schild_state(N=9, L=20.0, M=1.0, a=1.1)


def test_det_gtilde_unity():
    s = kerr_schild_state(N=9, L=20.0, M=1.0, a=0.5)
    dets = np.linalg.det(s.gtilde)
    assert np.allclose(dets, 1.0, atol=1e-10)
