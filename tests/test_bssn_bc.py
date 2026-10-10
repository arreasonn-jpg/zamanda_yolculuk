import numpy as np
from zt005.bssn_state import BSSNState, minkowski_state
from zt005.bssn_bc import sommerfeld_all


def test_sommerfeld_zero_when_no_change():
    s = minkowski_state(5)
    s_old = minkowski_state(5)
    out = sommerfeld_all(s, s_old, h=0.1, dt=0.01)
    for k, v in out.items():
        assert np.allclose(v, getattr(s, k)), f"{k} changed"


def test_sommerfeld_returns_all_keys():
    s = minkowski_state(5)
    out = sommerfeld_all(s, s, h=0.1, dt=0.01)
    for k in ("phi","gtilde","K","Atilde","Gamma","alpha","beta"):
        assert k in out


def test_sommerfeld_damps_wave_at_boundary():
    N = 5
    s = minkowski_state(N)
    s_old = minkowski_state(N)
    s.phi[-1, -1, -1] = 1.0
    out = sommerfeld_all(s, s_old, h=0.1, dt=0.05)
    # boundary value should be reduced compared to raw
    assert abs(out["phi"][-1,-1,-1]) < abs(s.phi[-1,-1,-1])
