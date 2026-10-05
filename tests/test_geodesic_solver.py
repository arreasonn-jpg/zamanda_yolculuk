"""Geodesic integrator validation against known answers (ZT-006.3).

The norm g(u,u) is an exact constant of motion along geodesics; drift is
the primary integrator diagnostic.  The Schwarzschild perihelion test
compares against the first post-Newtonian prediction
delta_phi = 6 pi M / (a (1 - e^2)).
"""
import numpy as np
import pytest

from zt005.benchmark_metrics import (GodelMetric, KerrMetric,
                                     MinkowskiMetric, MorrisThorneMetric,
                                     SchwarzschildMetric)
from zt005.geodesic_solver import (circular_orbit_initial_conditions_schwarzschild,
                                   integrate_geodesic)


def _normalize(metric, x, u, target=-1.0):
    n = float(u @ metric.g(x) @ u)
    assert n < 0, "test setup: initial 4-velocity must be timelike"
    return u * np.sqrt(target / n)


def test_minkowski_straight_line_exact():
    m = MinkowskiMetric()
    x0 = np.array([0.0, 0.0, 0.0, 0.0])
    u0 = np.array([2.0, 1.0, 0.5, -0.3])
    res = integrate_geodesic(m, x0, u0, tau_span=5.0, n_steps=20)
    assert np.allclose(res["x"][-1], x0 + 5.0 * u0, atol=1e-10)
    assert res["norm_drift_rel"] < 1e-12


def test_schwarzschild_circular_orbit_norm_and_radius():
    m = SchwarzschildMetric(M=1.0)
    x0, u0 = circular_orbit_initial_conditions_schwarzschild(1.0, 10.0)
    assert abs(float(u0 @ m.g(x0) @ u0) + 1.0) < 1e-12
    res = integrate_geodesic(m, x0, u0, tau_span=60.0, n_steps=4000,
                             analytic=True)
    assert res["norm_drift_rel"] < 1e-9
    # radius stays constant for a circular orbit
    assert np.max(np.abs(res["x"][:, 1] - 10.0)) < 1e-4


@pytest.mark.parametrize("metric, x0, u0_raw", [
    (KerrMetric(M=1.0, a=0.7), [0, 10.0, np.pi / 3, 0.0],
     [2.0, 0.1, 0.05, 0.05]),
    (GodelMetric(omega=1.0), [0, 0.5, 0.3, 0.2], [1.0, 0.2, 0.3, 0.1]),
    (MorrisThorneMetric(b0=1.0), [0, -3.0, np.pi / 2, 0.0],
     [1.0, 0.6, 0.0, 0.0]),
], ids=["kerr", "godel", "morris_thorne"])
def test_norm_conservation_numeric_connection(metric, x0, u0_raw):
    u0 = _normalize(metric, np.array(x0, float), np.array(u0_raw, float))
    res = integrate_geodesic(metric, x0, u0, tau_span=8.0, n_steps=1200,
                             numeric_step=1e-5)
    assert res["norm_drift_rel"] < 1e-5


def _turning_points(r, phi, which):
    idx = []
    for i in range(1, len(r) - 1):
        if which == "min" and r[i] < r[i - 1] and r[i] < r[i + 1]:
            idx.append(i)
        if which == "max" and r[i] > r[i - 1] and r[i] > r[i + 1]:
            idx.append(i)
    return idx


def test_schwarzschild_perihelion_precession_matches_1pn():
    """Bound eccentric orbit: measured advance per orbit vs 6piM/(a(1-e^2))."""
    M = 1.0
    m = SchwarzschildMetric(M=M)
    r_p = 40.0
    x0 = np.array([0.0, r_p, np.pi / 2, 0.0])
    omega_circ = np.sqrt(M / r_p ** 3)
    f = 1.0 - 2.0 * M / r_p
    # start exactly at perihelion: u^r = 0, tangential velocity slightly
    # above circular; normalization -f (u^t)^2 + r^2 (u^phi)^2 = -1
    k = 1.08 * omega_circ
    u_t = np.sqrt(1.0 / (f - r_p * r_p * k * k))
    u0 = np.array([u_t, 0.0, 0.0, k * u_t])
    assert abs(float(u0 @ m.g(x0) @ u0) + 1.0) < 1e-12

    res = integrate_geodesic(m, x0, u0, tau_span=2500.0, n_steps=50000,
                             analytic=True)
    r, phi = res["x"][:, 1], res["x"][:, 3]
    mins = _turning_points(r, phi, "min")
    assert len(mins) >= 1, "need the next perihelion"
    # the trajectory starts at a perihelion (u^r = 0, r increasing after)
    assert r[1] > r[0]
    dphi = float(phi[mins[0]] - phi[0])

    r_a = float(np.max(r[:mins[0] + 1]))
    a = 0.5 * (r_p + r_a)
    e = (r_a - r_p) / (r_a + r_p)
    predicted = 2.0 * np.pi + 6.0 * np.pi * M / (a * (1.0 - e * e))
    assert abs(dphi - predicted) / predicted < 0.05
