import numpy as np
import pytest
from zt005.sc_like_ic_kerr import sc_like_kerr_state
from zt005.sc_like_ic import sc_like_state


def test_a_zero_matches_nonrotating():
    s_k = sc_like_kerr_state(N=9, L=20.0, M=1.0, rc=2.0, a=0.0)
    s_n = sc_like_state(N=9, L=20.0, M=1.0, rc=2.0)
    assert np.allclose(s_k.phi, s_n.phi)
    assert np.allclose(s_k.beta, 0.0)


def test_a_nonzero_beta_nonzero():
    s = sc_like_kerr_state(N=9, L=20.0, M=1.0, rc=2.0, a=0.1)
    assert np.max(np.abs(s.beta)) > 1e-6


def test_beta_scales_with_a():
    s1 = sc_like_kerr_state(N=9, L=20.0, M=1.0, rc=2.0, a=0.1)
    s2 = sc_like_kerr_state(N=9, L=20.0, M=1.0, rc=2.0, a=0.2)
    r = np.max(np.abs(s2.beta)) / np.max(np.abs(s1.beta))
    assert abs(r - 2.0) < 0.01


def test_a_too_large_raises():
    with pytest.raises(ValueError):
        sc_like_kerr_state(N=9, L=20.0, M=1.0, rc=2.0, a=0.8)


def test_alpha_unchanged_by_rotation():
    s_k = sc_like_kerr_state(N=9, L=20.0, M=1.0, rc=2.0, a=0.3)
    s_n = sc_like_state(N=9, L=20.0, M=1.0, rc=2.0)
    assert np.allclose(s_k.alpha, s_n.alpha)
