"""Retarded Green-function solver tests (ZT-006.3 BLOCKER item).

1. Causality: a pulse emitted at the origin is observed at distance R
   only at t >= t0 + R/c, with the peak at the retarded time.
2. Static limit: a time-independent source reproduces the instantaneous
   integral.
"""
import numpy as np

from zt005.manufactured_sources import GaussianBump, conserved_source_fn
from zt005.retarded_gr import C, KAPPA, retarded_hbar, static_hbar


def test_retarded_pulse_arrives_at_retarded_time():
    t0, sigma_t = 0.0, 1e-8
    R = 30.0  # m
    w = 1.0   # m^3 single source cell

    def point_pulse(t, x):
        T = np.zeros((len(np.atleast_1d(t)), 4, 4))
        T[:, 0, 0] = np.exp(-(np.atleast_1d(t) - t0) ** 2 / (2 * sigma_t ** 2))
        return T

    obs = np.array([[R, 0.0, 0.0]])
    times = np.linspace(-4e-8, 1.4e-7, 361)
    h = retarded_hbar(point_pulse, obs, times,
                      np.array([[0.0, 0.0, 0.0]]), np.array([w]),
                      trace_reverse_source=False)
    signal = h[:, 0, 0, 0] / (KAPPA * w / R)   # normalized analytic envelope

    t_arrival = t0 + R / C
    peak_time = times[int(np.argmax(signal))]
    assert abs(peak_time - t_arrival) < (times[1] - times[0])
    # before the light front: negligible (Gaussian at 6 sigma ~ e^-18)
    before = np.max(np.abs(signal[times < t_arrival - 6 * sigma_t]))
    assert before < 1e-4
    # peak height: exact in the continuum; finite time sampling only
    assert abs(np.max(signal) - 1.0) < 1e-3


def test_static_limit_matches_instantaneous_integral():
    bump = GaussianBump(A=1.0, t0=0.0, sigma_t=1.0, sigma_x=1.0)
    src = GaussianBump(A=1.0, t0=0.0, sigma_t=1e30, sigma_x=1.0)  # frozen

    grid = np.linspace(-3.0, 3.0, 13)
    X, Y, Z = np.meshgrid(grid, grid, grid, indexing="ij")
    pts = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=1)
    vol = (grid[1] - grid[0]) ** 3
    weights = np.full(len(pts), vol)

    obs = np.array([[5.0, 0.0, 0.0], [0.0, 4.0, 0.0]])
    T_src = conserved_source_fn(src)(np.zeros(len(pts)), pts)

    h_ret = retarded_hbar(conserved_source_fn(src), obs, [0.0],
                          pts, weights)
    h_stat = static_hbar(T_src, obs, pts, weights)
    rel_err = float(np.max(np.abs(h_ret[0] - h_stat))
                    / np.max(np.abs(h_stat)))
    assert rel_err < 1e-9
