import numpy as np
from zt005.bssn_state_mp import minkowski_mp_state
from zt005.bssn_rk4_mp_vec_damp import rk4_step_mp_vec_damp

def test_minkowski_stationary():
    s = minkowski_mp_state(5)
    s2 = rk4_step_mp_vec_damp(s, h=0.1, dt=1e-5, width_cells=1, kappa=0.5)
    assert np.allclose(s2.phi, 0.0, atol=1e-10)
    assert np.allclose(s2.K, 0.0, atol=1e-10)

def test_puncture_finite():
    from zt005.bssn_ic_puncture import schwarzschild_puncture_state
    from zt005.bssn_state_mp import BSSNStateMP
    N = 7; L = 40.0; h = L/(N-1)
    s = BSSNStateMP.from_base(
        schwarzschild_puncture_state(N, L, M=1.0, eps=1.0))
    s2 = rk4_step_mp_vec_damp(s, h=h, dt=1e-4)
    assert np.all(np.isfinite(s2.phi))
    assert np.all(np.isfinite(s2.K))
