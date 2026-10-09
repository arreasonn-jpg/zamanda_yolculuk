"""Tests for exotic-matter sources and energy-gap report."""
import math
import numpy as np
import pytest

from zt005.exotic_matter import (
    casimir_parallel_plates,
    morris_thorne_shell,
    tipler_cylinder,
    compare_to_em_source,
)
from zt005.energy_gap_report import (
    build_report,
    energy_for_metric_perturbation,
    planck_scale_energy_J,
    sun_total_energy_J,
)


def test_casimir_negative_density():
    s = casimir_parallel_plates(separation_m=1e-6)
    assert s.rho(1e-6) < 0
    # d^-4 scaling
    s2 = casimir_parallel_plates(separation_m=2e-6)
    ratio = s2.rho(1e-6) / s.rho(1e-6)
    assert math.isclose(ratio, 1/16, rel_tol=1e-12)


def test_morris_thorne_negative_density():
    s = morris_thorne_shell(b0_m=1.0, throat_radius_m=1.0)
    assert s.rho(1.0) < 0
    # r^-4 scaling away from throat
    r1, r2 = 2.0, 4.0
    ratio = s.rho(r2) / s.rho(r1)
    assert math.isclose(ratio, (r1/r2)**4, rel_tol=1e-6)


def test_tipler_negative_density():
    s = tipler_cylinder(radius_m=1.0, omega_rad_s=1e6)
    assert s.rho(1.0) < 0
    # rho ~ -omega^2 R^2 r^2 / (8 pi G): quadratic in r
    r1, r2 = 1.0, 2.0
    ratio = s.rho(r2) / s.rho(r1)
    assert math.isclose(ratio, 4.0, rel_tol=1e-9)


def test_compare_to_em_source_positive_log():
    s = casimir_parallel_plates(separation_m=1e-6)
    d = compare_to_em_source(1e-30, s, r_probe=1e-6)
    assert d["log10_ratio"] > 0
    assert d["exotic_rho_J_per_m3"] < 0


def test_energy_for_perturbation_monotonic_in_V():
    E1 = energy_for_metric_perturbation(1e-3, 1e-3)
    E2 = energy_for_metric_perturbation(1e-3, 1.0)
    E3 = energy_for_metric_perturbation(1e-3, 1e3)
    assert E1 < E2 < E3


def test_report_has_all_fields():
    r = build_report()
    assert "h_gap_orders" in r
    assert r["h_gap_orders"]["gap_orders"] == pytest.approx(43.0, abs=0.1)
    assert "casimir_d1um" in r["exotic_vs_em"]
    assert "morris_thorne_b0_1m" in r["exotic_vs_em"]
    assert "tipler_1m_1e6rad_s" in r["exotic_vs_em"]


def test_planck_and_sun_scale_sane():
    assert 1e9 < planck_scale_energy_J() < 1e10
    assert 1e47 < sun_total_energy_J() < 1e48
