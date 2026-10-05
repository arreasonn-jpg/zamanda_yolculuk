
"""ZT-006.2.23 tests."""
import numpy as np
from zt005.run_zt006_2_23 import (
    extract_segments, biot_savart_B, compute_T00_physical,
)
from zt005.finite_wire_sources import build_finite_sources_grouped
from zt005.model import Backpack
from zt005.poisson_solver import solve_poisson_axisym


def test_extract_segments_nonempty():
    b = Backpack()
    s = build_finite_sources_grouped("G1", b, radius_m=0.005)
    mids, dls, I = extract_segments(s, I_wire=100.0)
    assert len(mids) > 100
    assert mids.shape == dls.shape
    assert len(I) == len(mids)
    assert np.all(I > 0)


def test_biot_savart_single_loop_on_axis():
    # One circular loop in xy-plane at z=0, radius 0.1, I=1 A
    n = 128
    phi = np.linspace(0, 2*np.pi, n, endpoint=False)
    pts = np.stack([0.1*np.cos(phi), 0.1*np.sin(phi), np.zeros(n)], axis=1)
    nxt = np.roll(pts, -1, axis=0)
    mids = 0.5*(pts+nxt)
    dls  = nxt - pts
    I    = np.ones(n) * 1.0
    B = biot_savart_B(mids, dls, I, np.array([0.0]), np.array([0.0]),
                      np.array([0.0]))
    # On-axis B = mu0 I / (2 R) = 4pi e-7 / 0.2 ~ 6.28e-6 T
    assert abs(B[0,2] - 6.283e-6) / 6.283e-6 < 0.02


def test_T00_physical_positive_for_all():
    r = np.linspace(0, 1, 20)
    z = np.linspace(-1.25, 1.25, 40)
    for g in ["G1","G2","G3","G4"]:
        T00 = compute_T00_physical(g, r, z, n_phi=4)
        assert np.max(T00) > 0, f"{g} gave zero T00"


def test_T00_physical_in_range():
    """Peak T00 should be ~1e-3..1e2 J/m^3 for 100 A / mm wires."""
    r = np.linspace(0, 1, 20)
    z = np.linspace(-1.25, 1.25, 40)
    T00 = compute_T00_physical("G1", r, z, n_phi=4)
    peak = np.max(T00)
    assert 1e-6 < peak < 1e4, f"Peak T00 = {peak} not in plausible range"


def test_poisson_from_physical_source_finite():
    r = np.linspace(0, 1, 20)
    z = np.linspace(-1.25, 1.25, 40)
    T00 = compute_T00_physical("G1", r, z, n_phi=4)
    h = solve_poisson_axisym(T00, r, z)
    assert np.all(np.isfinite(h))
    assert np.max(np.abs(h)) > 0
