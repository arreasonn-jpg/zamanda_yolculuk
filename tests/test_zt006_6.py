import pytest
pytestmark = pytest.mark.slow

"""ZT-006.6 tests - device gravitomagnetic field."""
import numpy as np
from zt005.run_zt006_6 import (
    compute_A_phi_wire, compute_T_0phi_peak, run_geometry,
)
from zt005.finite_wire_sources import build_finite_sources_grouped
from zt005.run_zt006_2_23 import extract_segments
from zt005.model import Backpack


def test_A_phi_single_loop_positive():
    """On-axis A_phi = 0; off-axis A_phi > 0 for z-directed current loop."""
    n = 128
    phi = np.linspace(0, 2*np.pi, n, endpoint=False)
    pts = np.stack([0.1*np.cos(phi), 0.1*np.sin(phi), np.zeros(n)], axis=1)
    nxt = np.roll(pts, -1, axis=0)
    mids = 0.5*(pts+nxt)
    dls  = nxt - pts
    I    = np.ones(n) * 1.0
    A_at_r = compute_A_phi_wire(mids, dls, I, 0.1, 0.0)
    A_far = compute_A_phi_wire(mids, dls, I, 1.0, 0.0)
    # A_phi should be positive and larger near the loop
    assert A_at_r > 0
    assert A_far > 0
    assert A_at_r > A_far


def test_T_0phi_zero_for_axisym_solenoid():
    """For axisymmetric solenoid, T_0phi = 0 analytically."""
    rv = np.linspace(0.0, 1.0, 10)
    zv = np.linspace(-1.25, 1.25, 20)
    T = compute_T_0phi_peak("G1", rv, zv, n_phi=4)
    assert T.shape == (10, 20)
    assert np.allclose(T, 0.0, atol=1e-100)


def test_T_0phi_zero_for_all_geometries():
    """All current geometries give T_0phi = 0 due to axisym + E_phi-only."""
    rv = np.linspace(0.0, 1.0, 8)
    zv = np.linspace(-1.0, 1.0, 16)
    for geom in ["G1", "G2", "G3", "G4"]:
        T = compute_T_0phi_peak(geom, rv, zv, n_phi=4)
        assert np.allclose(T, 0.0, atol=1e-100), f"{geom} not zero"


def test_run_geometry_returns_finite():
    rv, zv, T, h = run_geometry("G1", Nr=20, Nz=40, n_phi=4)
    assert np.all(np.isfinite(T))
    assert np.all(np.isfinite(h))
