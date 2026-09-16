"""Validation test suite for 5 benchmark GR spacetimes."""

import pytest
from zt005.benchmark_suite import (
    run_minkowski_benchmark,
    run_schwarzschild_benchmark,
    run_kerr_benchmark,
    run_godel_benchmark,
    run_tipler_benchmark,
    run_full_benchmark_suite,
)


def test_minkowski():
    res = run_minkowski_benchmark()
    assert res["pass"] is True
    assert res["metric_error"] == 0.0
    assert res["christoffel_error"] == 0.0


def test_schwarzschild():
    res = run_schwarzschild_benchmark(M=1.0)
    assert res["pass"] is True
    assert res["metric_error"] == 0.0
    assert res["curvature_error"]["ricci_vacuum_residual"] < 1e-5


def test_kerr():
    res = run_kerr_benchmark(M=1.0, a=0.7)
    assert res["pass"] is True
    assert res["metric_symmetry_error"] == 0.0
    assert res["curvature_error"]["ricci_vacuum_residual"] < 1e-5


def test_godel():
    res = run_godel_benchmark(omega=1.0)
    assert res["pass"] is True
    assert res["ctc_behavior_correct"] is True


def test_tipler():
    res = run_tipler_benchmark(a=1.0)
    assert res["pass"] is True
    assert res["ctc_behavior_correct"] is True


def test_full_benchmark_suite():
    res = run_full_benchmark_suite()
    assert res["overall_pass"] is True
