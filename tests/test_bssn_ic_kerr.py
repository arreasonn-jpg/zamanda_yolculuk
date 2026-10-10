import numpy as np
import pytest
from zt005.bssn_ic_kerr import kerr_slow_rotation_state
from zt005.bssn_ic_puncture import schwarzschild_puncture_state


def test_a_zero_matches_schwarzschild():
    s_k = kerr_slow_rotation_state(N=9, L=20.0, M=1.0, a=0.0)
    s_s = schwarzschild_puncture_state(N=9, L=20.0, M=1.0)
    assert np.allclose(s_k.phi, s_s.phi)
    assert np.allclose(s_k.alpha, s_s.alpha)
    assert np.allclose(s_k.beta, 0.0)


def test_a_nonzero_beta_nonzero():
    s = kerr_slow_rotation_state(N=9, L=20.0, M=1.0, a=0.1)
    assert np.max(np.abs(s.beta)) > 1e-6


def test_beta_scales_with_a():
    s1 = kerr_slow_rotation_state(N=9, L=20.0, M=1.0, a=0.1)
    s2 = kerr_slow_rotation_state(N=9, L=20.0, M=1.0, a=0.2)
    r = np.max(np.abs(s2.beta)) / np.max(np.abs(s1.beta))
    assert abs(r - 2.0) < 0.01


def test_a_too_large_raises():
    with pytest.raises(ValueError):
        kerr_slow_rotation_state(N=9, L=20.0, M=1.0, a=0.8)
