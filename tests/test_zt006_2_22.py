
"""ZT-006.2.22 unit tests."""
import numpy as np
from zt005.run_zt006_2_22 import (
    _flatten_sources,
    build_T00_axisym_smoothed,
)
from zt005.finite_wire_sources import build_finite_sources_grouped
from zt005.poisson_solver import solve_poisson_axisym
from zt005.model import Backpack


def test_flatten_sources_returns_nonempty_points():
    b = Backpack()
    sources = build_finite_sources_grouped("G1", b, radius_m=0.005)
    pts, currents = _flatten_sources(sources)
    assert pts.ndim == 2
    assert pts.shape[1] == 3
    assert len(pts) > 100, f"Expected >100 points, got {len(pts)}"
    assert len(currents) == len(pts)
    assert np.all(currents > 0)


def test_T00_nonzero_for_all_geometries():
    b = Backpack()
    r = np.linspace(0, 1, 30)
    z = np.linspace(-1.25, 1.25, 60)
    for geom in ["G1", "G2", "G3", "G4"]:
        T00 = build_T00_axisym_smoothed(geom, r, z, backpack=b)
        assert np.max(T00) > 0, f"{geom} produced zero T00"


def test_T00_shapes_match_grid():
    r = np.linspace(0, 1, 25)
    z = np.linspace(-1.25, 1.25, 50)
    T00 = build_T00_axisym_smoothed("G1", r, z)
    assert T00.shape == (25, 50)


def test_poisson_from_real_source_is_finite():
    r = np.linspace(0, 1, 30)
    z = np.linspace(-1.25, 1.25, 60)
    T00 = build_T00_axisym_smoothed("G1", r, z)
    h = solve_poisson_axisym(T00, r, z)
    assert np.all(np.isfinite(h))
