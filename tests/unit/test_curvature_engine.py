"""Unit tests for full curvature and Einstein tensor pipeline."""

import numpy as np
import pytest
from zt005.curvature_engine import (
    compute_metric_inverse,
    compute_christoffel,
    compute_riemann,
    compute_curvature_chain,
    compute_einstein_residual,
)
from zt005.benchmark_metrics import MinkowskiMetric, SchwarzschildMetric


def test_metric_inverse():
    eta = np.diag([-1.0, 1.0, 1.0, 1.0])
    inv = compute_metric_inverse(eta)
    assert np.allclose(inv, eta)


def test_minkowski_flatness():
    m = MinkowskiMetric()
    probe = np.array([0.0, 1.0, 2.0, 3.0])
    chain = compute_curvature_chain(m.g, probe)
    assert np.allclose(chain["gamma"], 0.0, atol=1e-12)
    assert np.allclose(chain["riemann"], 0.0, atol=1e-12)
    assert np.allclose(chain["ricci"], 0.0, atol=1e-12)
    assert abs(chain["ricci_scalar"]) < 1e-12
    assert np.allclose(chain["einstein"], 0.0, atol=1e-12)


def test_schwarzschild_vacuum():
    m = SchwarzschildMetric(M=1.0)
    probe = np.array([0.0, 6.0, np.pi / 2.0, 0.0])
    chain = compute_curvature_chain(m.g, probe, h=1e-4)
    # Riemann is non-zero (tidal gravity)
    assert np.linalg.norm(chain["riemann"]) > 1e-4
    # Ricci and Einstein must vanish in vacuum
    assert np.allclose(chain["ricci"], 0.0, atol=1e-5)
    assert np.allclose(chain["einstein"], 0.0, atol=1e-5)


def test_einstein_residual_solver_classification():
    m = MinkowskiMetric()
    probe = np.array([0.0, 0.0, 0.0, 0.0])
    res = compute_einstein_residual(m.g, probe, np.zeros((4, 4)))
    assert res["is_full_gr_solution"] is True
    assert res["classification"] == "linearized metric perturbation solver"
