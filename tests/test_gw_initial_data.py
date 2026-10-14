import numpy as np
import pytest
from zt005.gw_initial_data import gw_packet

def test_packet_shapes():
    out = gw_packet(N=11, L=20.0, A_amp=1e-4, sigma=1.0, z0=0.0)
    for k in ("gtilde", "phi", "Atilde", "K", "h_ij"):
        assert out[k].shape[:3] == (11, 11, 11)

def test_plus_structure():
    out = gw_packet(N=11, L=20.0, A_amp=1e-4, sigma=1.0, z0=0.0,
                    polarization="plus")
    h = out["h_ij"]
    mid = 5
    assert h[mid, mid, mid, 0, 0] > 0
    assert h[mid, mid, mid, 1, 1] < 0
    assert abs(h[mid, mid, mid, 0, 1]) < 1e-15

def test_cross_structure():
    out = gw_packet(N=11, L=20.0, A_amp=1e-4, sigma=1.0, z0=0.0,
                    polarization="cross")
    h = out["h_ij"]
    mid = 5
    assert abs(h[mid, mid, mid, 0, 1]) > 1e-5

def test_det_gtilde_unity():
    out = gw_packet(N=11, L=20.0, A_amp=1e-4, sigma=1.0, z0=0.0)
    dets = np.linalg.det(out["gtilde"])
    assert np.allclose(dets, 1.0, atol=1e-10)

def test_amplitude_scaling():
    o1 = gw_packet(N=11, L=20.0, A_amp=1e-4, sigma=1.0, z0=0.0)
    o2 = gw_packet(N=11, L=20.0, A_amp=2e-4, sigma=1.0, z0=0.0)
    assert np.max(np.abs(o2["h_ij"])) == pytest.approx(
        2.0 * np.max(np.abs(o1["h_ij"])), rel=1e-10)
