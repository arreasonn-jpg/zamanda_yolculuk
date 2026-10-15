import numpy as np
import pytest
from zt005.teleparallel import (
    tetrad_flat, torsion_scalar_flat,
    f_T_gravity, f_T_GR_limit, ptum_qft_window_in_fT,
)

def test_flat_tetrad_identity():
    e = tetrad_flat()
    assert np.allclose(e, np.eye(4))

def test_flat_torsion_zero():
    assert torsion_scalar_flat() == 0.0

def test_fT_GR_limit():
    """alpha=0 -> f(T) = T (TEGR)."""
    assert f_T_gravity(1.0, alpha=0.0) == 1.0

def test_fT_modification():
    """alpha>0 -> deviation from TEGR."""
    assert f_T_gravity(1.0, alpha=1.0) == 2.0

def test_GR_equivalence_marker():
    g = f_T_GR_limit()
    assert g["equivalent_to_GR"]

def test_ptum_in_fT_invariant():
    """PTUM QFT window unchanged."""
    w1 = ptum_qft_window_in_fT(alpha=0.0)
    w2 = ptum_qft_window_in_fT(alpha=0.1)
    assert w1["A_min"] == w2["A_min"] == 1.35
    assert w1["A_max"] == w2["A_max"] == 1.87
