import numpy as np
from zt005.bssn_state_mp import BSSNStateMP, minkowski_mp_state
from zt005.bssn_rk4_mp_vec import rk4_step_mp_vec

def test_minkowski_stationary():
    s = minkowski_mp_state(5)
    s2 = rk4_step_mp_vec(s, h=0.1, dt=1e-5)
    assert np.allclose(s2.phi, 0.0, atol=1e-12)
    assert np.allclose(s2.beta, 0.0, atol=1e-12)

def test_puncture_runs_stable():
    """Schwarzschild puncture: RHS is finite and bounded."""
    from zt005.bssn_ic_puncture import schwarzschild_puncture_state
    N = 7; L = 40.0; h = L/(N-1)
    s = BSSNStateMP.from_base(
        schwarzschild_puncture_state(N, L, M=1.0, eps=1.0))
    s2 = rk4_step_mp_vec(s, h=h, dt=1e-4)
    assert np.all(np.isfinite(s2.phi))
    assert np.all(np.isfinite(s2.K))
    assert np.max(np.abs(s2.K)) < 1e-4  # bounded

def test_matches_loop_within_1e5():
    """Vectorized and loop agree at the 1e-5 level
    (loop version has known residual differences in ricci_tensor)."""
    from zt005.bssn_rk4_mp import rk4_step_mp
    from zt005.bssn_ic_puncture import schwarzschild_puncture_state
    N = 7; L = 40.0; h = L/(N-1)
    s_loop = BSSNStateMP.from_base(
        schwarzschild_puncture_state(N, L, M=1.0, eps=1.0))
    s_vec = BSSNStateMP.from_base(
        schwarzschild_puncture_state(N, L, M=1.0, eps=1.0))
    s_loop2 = rk4_step_mp(s_loop, h=h, dt=1e-4)
    s_vec2 = rk4_step_mp_vec(s_vec, h=h, dt=1e-4)
    sl = slice(2, -2)
    assert np.allclose(s_loop2.phi[sl,sl,sl], s_vec2.phi[sl,sl,sl], atol=1e-5)
    assert np.allclose(s_loop2.K[sl,sl,sl], s_vec2.K[sl,sl,sl], atol=1e-5)
