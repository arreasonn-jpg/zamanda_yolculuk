"""Lorenz-gauge consistency tests (ZT-006.3 critical item 4).

The exact retarded solution of a CONSERVED source satisfies
partial_mu hbar^{mu nu} = 0 identically.  Numerically: build the retarded
solution with retarded_gr.retarded_hbar from the manufactured conserved
source, then differentiate it on a spacetime grid.  The residual must be
small relative to hbar / (c sigma_t), the natural gradient scale.
"""
import numpy as np
import pytest

from zt005.gauge_checks import lorenz_residual
from zt005.manufactured_sources import GaussianBump, conserved_source_fn
from zt005.retarded_gr import C, retarded_hbar


def _retarded_hbar_grid():
    bump = GaussianBump(A=1.0, t0=0.0, x0=(0.0, 0.0, 0.0),
                        sigma_t=1e-7, sigma_x=1.0)

    # source quadrature grid covering the bump support.  The spacing must
    # resolve the source well because quadrature noise is amplified by the
    # finite-difference step of the Lorenz residual.
    s = np.linspace(-3.0, 3.0, 21)
    X, Y, Z = np.meshgrid(s, s, s, indexing="ij")
    src_pts = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=1)
    vol = (s[1] - s[0]) ** 3
    weights = np.full(len(src_pts), vol)

    # observation grid (interior points are used for differences).
    # spacing dx = 0.3 m keeps the FD truncation error of the residual
    # itself small: O(dx^2) * d^3 hbar.
    o = np.linspace(-0.6, 0.6, 7)
    OX, OY, OZ = np.meshgrid(o, o, o, indexing="ij")
    obs_pts = np.stack([OX.ravel(), OY.ravel(), OZ.ravel()], axis=1)

    times = np.linspace(-6e-8, 6e-8, 13)
    h = retarded_hbar(conserved_source_fn(bump), obs_pts, times,
                      src_pts, weights)
    nt, no = len(times), len(o)
    grid = h.reshape(nt, no, no, no, 4, 4)
    return grid, times[1] - times[0], o[1] - o[0], bump


@pytest.fixture(scope="module")
def residual_report():
    grid, dt, dx, bump = _retarded_hbar_grid()
    rep = lorenz_residual(grid, dt, dx)
    # natural gradient scale: variation happens over ~c*sigma_t or sigma_x
    L = min(C * bump.st, bump.sx)
    rep["scaled_residual"] = rep["max_abs"] * L / rep["scale"]
    return rep


def test_lorenz_residual_small_for_conserved_source(residual_report):
    rep = residual_report
    assert rep["scaled_residual"] < 0.15, rep["scaled_residual"]


def test_lorenz_residual_structurally_nonzero(residual_report):
    # sanity: the solution itself is not identically zero
    assert residual_report["scale"] > 0.0
