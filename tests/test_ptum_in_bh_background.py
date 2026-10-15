import numpy as np
import pytest
from zt005.ptum_in_bh_background import (
    schwarzschild_f, ptum_static_in_bh, ptum_qft_threshold_in_bh,
)

def test_schwarzschild_far_field():
    """r >> M -> f -> 1."""
    assert schwarzschild_f(100.0, M=1.0) == pytest.approx(1.0, abs=3e-2)

def test_schwarzschild_horizon_zero():
    """At r = M/2 (horizon), f = 0."""
    f = schwarzschild_f(0.5, M=1.0)
    assert f < 1e-6

def test_ptum_activation_near_peak():
    """Near Gaussian peak (r=1): PTUM activated."""
    r = ptum_static_in_bh(M=1.0, r_probe=1.0, A_fuel=3.0)
    assert r["meff2"] < 0
    assert r["psi_eq"] > 0
    assert r["V_min"] < 0

def test_ptum_deactivated_far_from_peak():
    """Far from Gaussian peak (r=8): Phi < Phi_c, no activation."""
    r = ptum_static_in_bh(M=1.0, r_probe=8.0, A_fuel=3.0)
    assert r["meff2"] > 0
    assert r["psi_eq"] == 0.0
    assert r["V_min"] == 0.0

def test_qft_threshold_invariant():
    """QFT window independent of M, r_probe."""
    w1 = ptum_qft_threshold_in_bh(M=1.0, r_probe=10.0)
    w2 = ptum_qft_threshold_in_bh(M=5.0, r_probe=2.0)
    assert w1["A_min"] == w2["A_min"] == 1.35
    assert w1["A_max"] == w2["A_max"] == 1.87
    assert w1["invariant"]
