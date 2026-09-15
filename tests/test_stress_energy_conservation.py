"""partial_mu T^{mu nu} = 0 diagnostic tests (ZT-006.3 critical item 3).

The conserved manufactured source is built from the linearized Einstein
tensor of a pure-trace bump, which is identically divergence-free; the
broken source is a control that any honest diagnostic must flag.
"""
import numpy as np

from zt005.manufactured_sources import (GaussianBump, broken_source_fn,
                                        conserved_source_fn)
from zt005.stress_energy_conservation import (divergence_residual,
                                              sample_source_on_grid)


def _grid(bump):
    st, sx = bump.st, bump.sx
    t = np.linspace(bump.t0 - 3 * st, bump.t0 + 3 * st, 25)
    x = np.linspace(-3 * sx, 3 * sx, 25)
    return t, x, x, x


def test_conserved_source_residual_small():
    bump = GaussianBump(A=1.0, sigma_t=1e-7, sigma_x=1.0)
    t, x, y, z = _grid(bump)
    T = sample_source_on_grid(conserved_source_fn(bump), t, x, y, z)
    dt = t[1] - t[0]
    dx = x[1] - x[0]
    res = divergence_residual(T, dt, dx)
    # second-order finite differences on a smooth, exactly-conserved source:
    # residual must be far below the source scale (FD truncation only)
    assert res["max_rel"] < 2e-2, res["max_rel"]


def test_broken_source_residual_large():
    bump = GaussianBump(A=1.0, sigma_t=1e-7, sigma_x=1.0)
    t, x, y, z = _grid(bump)
    T = sample_source_on_grid(broken_source_fn(bump), t, x, y, z)
    res = divergence_residual(T, t[1] - t[0], x[1] - x[0])
    # violation at the spatial-gradient scale: order-unity relative signal
    assert res["max_rel"] > 0.2, res["max_rel"]


def test_residual_converges_with_resolution():
    bump = GaussianBump(A=1.0, sigma_t=1e-7, sigma_x=1.0)

    def maxrel(n):
        t = np.linspace(bump.t0 - 3 * bump.st, bump.t0 + 3 * bump.st, n)
        x = np.linspace(-3 * bump.sx, 3 * bump.sx, n)
        T = sample_source_on_grid(conserved_source_fn(bump), t, x, x, x)
        return divergence_residual(T, t[1] - t[0], x[1] - x[0])["max_rel"]

    r1, r2 = maxrel(21), maxrel(33)
    assert r2 < r1 / 2.5   # ~2nd order convergence
