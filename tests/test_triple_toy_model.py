import numpy as np
import pytest
from zt005.triple_toy_model import (
    sc_core_metric, ptum_local_activation, triple_state_at_r,
    qft_window_triple, triple_coupling_check,
)

def test_sc_core_regular_at_origin():
    """SC M_eff(0) = 0 (regular core)."""
    assert sc_core_metric(0.0, M=1.0, rc=2.0) == pytest.approx(0.0, abs=1e-10)

def test_sc_core_approaches_M_at_large_r():
    """Far from core: M_eff -> M."""
    assert sc_core_metric(20.0, M=1.0, rc=2.0) == pytest.approx(1.0, abs=1e-3)

def test_ptum_activation_above_threshold():
    assert ptum_local_activation(Phi=3.0)  # meff2 < 0
    assert not ptum_local_activation(Phi=0.5)  # meff2 > 0

def test_triple_state_at_r_shapes():
    s = triple_state_at_r(r=5.0, M=1.0, rc=2.0, A_fuel=3.0)
    for k in ("M_eff", "f_bh", "Phi", "psi_eq", "V_min", "activated"):
        assert k in s

def test_qft_window_invariant_triple():
    w = qft_window_triple()
    assert w["A_min"] == 1.35
    assert w["A_max"] == 1.87
    assert w["SC_shift"] == 0.0
    assert w["BH_shift"] == 0.0
    assert w["invariant"]

def test_coupling_check():
    """SC + BH present, PTUM activates locally (matter sector)."""
    c = triple_coupling_check(M=1.0, rc=2.0, r_probe=5.0)
    assert c["SC_M_eff"] > 0
    assert c["BH_lapse"] > 0
    # At r=5, Phi ~ 3*exp(-25/32) ~ 1.37 > 1 -> activated
    assert c["PTUM_activated"]
