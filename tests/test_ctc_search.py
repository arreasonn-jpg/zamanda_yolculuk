"""Known-answer CTC benchmark tests (ZT-006.3, PHASE 6 groundwork).

Positive answers: Gödel (r > asinh(1)) and the idealized Tipler cylinder
(a r > 1).  Negative answers: Minkowski, Schwarzschild exterior,
Morris-Thorne throat.  These tests validate the detector only; they make
no statement about the device model.
"""
import numpy as np

from zt005.benchmark_metrics import GodelMetric, TiplerCylinderMetric
from zt005.ctc_search import (check_phi_circle_ctc, godel_ctc_benchmark,
                              negative_benchmarks, tipler_ctc_benchmark)


def test_godel_benchmark_agrees_with_analytic():
    rep = godel_ctc_benchmark(omega=1.0, n_points=129)
    assert rep["agrees_with_analytic"]
    assert rep["detected_outside"] is True
    assert rep["detected_inside"] is False


def test_tipler_benchmark_agrees_with_analytic():
    rep = tipler_ctc_benchmark(a=1.0, n_points=129)
    assert rep["agrees_with_analytic"]
    assert rep["detected_outside"] is True
    assert rep["detected_inside"] is False


def test_negative_benchmarks_find_no_ctc():
    rep = negative_benchmarks(n_points=65)
    assert rep["all_clean"]
    assert rep["minkowski"]["found"] is False
    assert rep["schwarzschild_r6M"]["found"] is False
    assert rep["morris_thorne_throat"]["found"] is False


def test_godel_circle_proper_time_positive_inside_ctc_region():
    m = GodelMetric(omega=1.0)
    rep = check_phi_circle_ctc(m, 1.5 * m.ctc_radius(), coord_layout="trpz")
    assert rep["is_ctc_candidate"]
    assert rep["proper_time"] > 0.0
