import numpy as np
from zt005.bssn_rk4_grid import rk4_step_grid


def _minkowski_state(N):
    return {
        "phi": np.zeros((N,N,N)),
        "gt": np.broadcast_to(np.eye(3)[None,None,None,:,:], (N,N,N,3,3)).copy(),
        "K": np.zeros((N,N,N)),
        "At": np.zeros((N,N,N,3,3)),
    }


def test_rk4_step_minkowski_stationary():
    s = _minkowski_state(5)
    out = rk4_step_grid(s, h=0.1, dt=1e-3)
    assert np.allclose(out["phi"], 0.0, atol=1e-12)
    assert np.allclose(out["K"], 0.0, atol=1e-12)
    assert np.allclose(out["gt"], np.eye(3))


def test_rk4_step_returns_all_keys():
    s = _minkowski_state(5)
    out = rk4_step_grid(s, h=0.1, dt=1e-3)
    for k in ("phi", "gt", "K", "At"):
        assert k in out
        assert out[k].shape == s[k].shape


def test_rk4_step_phi_perturbation_finite():
    s = _minkowski_state(5)
    s["phi"][2,2,2] = 1e-5
    out = rk4_step_grid(s, h=0.1, dt=1e-4)
    assert np.all(np.isfinite(out["phi"]))
    assert np.max(np.abs(out["phi"])) < 1e-4
