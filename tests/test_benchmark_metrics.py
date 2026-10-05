"""Benchmark metric sanity tests (ZT-006.3)."""
import numpy as np
import pytest

from zt005.benchmark_metrics import (
    GodelMetric, KerrMetric, MinkowskiMetric, MorrisThorneMetric,
    SchwarzschildMetric, TiplerCylinderMetric, all_benchmarks)


def test_all_metrics_symmetric_and_lorentzian():
    probe = {
        "minkowski": np.array([0.0, 1.0, 1.0, 1.0]),
        "schwarzschild": np.array([0.0, 6.0, np.pi / 2, 0.0]),
        "kerr": np.array([0.0, 10.0, np.pi / 2, 0.0]),
        "godel": np.array([0.0, 1.2, 0.0, 0.0]),
        "tipler_cylinder": np.array([0.0, 1.5, 0.0, 0.0]),
        "morris_thorne": np.array([0.0, 0.0, np.pi / 2, 0.0]),
    }
    for m in all_benchmarks():
        g = m.g(probe[m.name])
        assert np.allclose(g, g.T), m.name
        ev = np.linalg.eigvalsh(g)
        assert np.sum(ev < 0) == 1 and np.sum(ev > 0) == 3, m.name


def test_minkowski_is_flat():
    m = MinkowskiMetric()
    assert np.allclose(m.g([0, 0, 0, 0]), np.diag([-1, 1, 1, 1]))
    assert np.allclose(m.christoffel_numeric([0, 0, 0, 0]), 0.0, atol=1e-10)


def test_schwarzschild_gtt_grr_product():
    m = SchwarzschildMetric(M=1.0)
    for r in (3.0, 6.0, 20.0):
        g = m.g([0.0, r, np.pi / 2, 0.0])
        assert abs(g[0, 0] * g[1, 1] + 1.0) < 1e-12


def test_schwarzschild_numeric_vs_analytic_christoffel():
    m = SchwarzschildMetric(M=1.0)
    x = np.array([0.0, 8.0, np.pi / 3, 0.7])
    ga = m.christoffel_analytic(x)
    gn = m.christoffel_numeric(x, h=1e-4)
    assert np.allclose(ga, gn, atol=1e-5)


def test_kerr_reduces_to_schwarzschild_at_a0():
    k = KerrMetric(M=1.0, a=1e-6)
    s = SchwarzschildMetric(M=1.0)
    x = [0.0, 9.0, np.pi / 3, 1.0]
    assert np.allclose(k.g(x), s.g(x), atol=1e-4)


def test_godel_ctc_radius_sign_change():
    m = GodelMetric(omega=1.0)
    r_G = m.ctc_radius()
    assert abs(r_G - np.arcsinh(1.0)) < 1e-14
    assert m.g([0, 0.5 * r_G, 0, 0])[2, 2] > 0      # spacelike phi-circles
    assert m.g([0, 1.5 * r_G, 0, 0])[2, 2] < 0      # timelike phi-circles
    # still Lorentzian even in the CTC region
    ev = np.linalg.eigvalsh(m.g([0, 1.5 * r_G, 0, 0]))
    assert np.sum(ev < 0) == 1


def test_tipler_ctc_radius_sign_change():
    m = TiplerCylinderMetric(a=1.0)
    assert m.g([0, 0.5, 0, 0])[2, 2] > 0
    assert m.g([0, 1.5, 0, 0])[2, 2] < 0
    # determinant of t-phi block is -r^2 (Lorentzian everywhere)
    g = m.g([0, 1.5, 0, 0])
    det_block = g[0, 0] * g[2, 2] - g[0, 2] ** 2
    assert det_block < 0


def test_morris_thorne_throat_radius():
    m = MorrisThorneMetric(b0=2.0)
    g = m.g([0.0, 0.0, np.pi / 2, 0.0])
    assert abs(g[2, 2] - 4.0) < 1e-12   # areal radius^2 = b0^2 at throat
    assert g[0, 0] == -1.0               # zero redshift function
