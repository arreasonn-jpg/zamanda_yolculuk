"""Causal classification tests — project convention: ds^2<0 timelike."""
import numpy as np

from zt005.benchmark_metrics import GodelMetric, MinkowskiMetric
from zt005.causal_structure import (NULL, SPACELIKE, TIMELIKE,
                                    classify_curve, classify_interval)


def test_interval_classification():
    assert classify_interval(-1.0) == TIMELIKE
    assert classify_interval(1.0) == SPACELIKE
    assert classify_interval(1e-15) == NULL


def test_minkowski_segment_classes():
    m = MinkowskiMetric()
    # timelike: motion slower than c (c=1 units)
    pts = np.array([[0, 0, 0, 0], [1, 0.5, 0, 0]], float)
    assert classify_curve(m, pts)["overall"] == TIMELIKE
    # null: 45 degrees
    pts = np.array([[0, 0, 0, 0], [1, 1.0, 0, 0]], float)
    assert classify_curve(m, pts)["overall"] == NULL
    # spacelike
    pts = np.array([[0, 0, 0, 0], [0.2, 1.0, 0, 0]], float)
    assert classify_curve(m, pts)["overall"] == SPACELIKE


def test_godel_phi_circle_classification_matches_analytic():
    m = GodelMetric(omega=1.0)
    r_G = m.ctc_radius()
    for r, expected in [(0.5 * r_G, SPACELIKE), (1.5 * r_G, TIMELIKE)]:
        phi = np.linspace(0, 2 * np.pi, 129)
        pts = np.zeros((129, 4))
        pts[:, 1] = r
        pts[:, 2] = phi
        rep = classify_curve(m, pts, periodic_coord=2)
        assert rep["overall"] == expected, (r, expected)
        assert rep["closed_covering_space"]
