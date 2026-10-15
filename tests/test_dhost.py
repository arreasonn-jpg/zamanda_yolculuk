import numpy as np
import pytest
from zt005.dhost import (
    dhost_parameters, beyond_horndeski_terms,
    gw_speed_ratio, gw170817_constraint,
    ptum_qft_window_in_dhost,
)

def test_dhost_GR_limit():
    """A1=A2=A3=0 -> GR."""
    p = dhost_parameters()
    assert p["A1"] == 0.0 and p["A2"] == 0.0 and p["A3"] == 0.0
    assert p["c_GW_deviation"] == 0.0

def test_gw_speed_GR_limit():
    """A1=0 -> c_GW = c."""
    assert gw_speed_ratio(0.0) == 1.0

def test_gw_speed_deviation():
    """A1 != 0 -> c_GW != c."""
    assert gw_speed_ratio(0.01) > 1.0
    assert gw_speed_ratio(-0.01) < 1.0

def test_gw170817_bound_tight():
    """GW170817 bound is extremely tight."""
    c = gw170817_constraint()
    assert c["bound"] == 1e-15

def test_beyond_horndeski_zero_at_GR():
    """d2phi=0 -> beyond terms zero."""
    d2phi = np.zeros((3, 3))
    assert beyond_horndeski_terms(1.0, 0.5, d2phi) == 0.0

def test_ptum_qft_invariant_in_dhost():
    """PTUM QFT window unchanged by DHOST."""
    w1 = ptum_qft_window_in_dhost(A1=0.0)
    w2 = ptum_qft_window_in_dhost(A1=0.1)
    w3 = ptum_qft_window_in_dhost(A1=-0.1)
    assert w1["A_min"] == w2["A_min"] == w3["A_min"] == 1.35
    assert w1["A_max"] == w2["A_max"] == w3["A_max"] == 1.87
