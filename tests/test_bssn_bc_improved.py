import numpy as np
import pytest
from zt005.bssn_bc_improved import (
    sponge_profile, sponge_3d, damp_K, multi_layer_sommerfeld,
)

def test_sponge_profile_interior_zero():
    p = sponge_profile(21, width_cells=3)
    assert p[10] == 0.0
    assert p[8] == 0.0

def test_sponge_profile_edges_one():
    p = sponge_profile(21, width_cells=3)
    assert p[0] == 1.0
    assert p[-1] == 1.0

def test_sponge_profile_smooth():
    p = sponge_profile(21, width_cells=5)
    # Monotonic from edge to interior
    for i in range(4):
        assert p[i] > p[i+1]

def test_sponge_3d_shape():
    P = sponge_3d(11, width_cells=2)
    assert P.shape == (11, 11, 11)

def test_sponge_3d_corners_one():
    P = sponge_3d(11, width_cells=2)
    assert P[0,0,0] == 1.0
    assert P[5,5,5] == 0.0

def test_damp_K_zero_when_sponge_zero():
    K = np.ones((5,5,5))
    H = np.ones((5,5,5))
    sponge = np.zeros((5,5,5))
    K2 = damp_K(K, H, sponge, kappa=1.0)
    assert np.allclose(K2, K)

def test_damp_K_reduces_near_edge():
    K = np.ones((5,5,5))
    H = np.ones((5,5,5))
    sponge = np.ones((5,5,5))
    K2 = damp_K(K, H, sponge, kappa=0.5)
    assert np.allclose(K2, 0.5)

def test_sommerfeld_no_change_for_static():
    f = np.ones((7,7,7))
    f_old = np.ones((7,7,7))
    out = multi_layer_sommerfeld(f, f_old, h=0.1, dt=0.01, width_cells=2)
    assert np.allclose(out, f)
