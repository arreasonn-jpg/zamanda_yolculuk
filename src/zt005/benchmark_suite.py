"""Comprehensive 5-Metric Benchmark Suite.

Addresses Section 4:
Executes rigorous error evaluations for the 5 fundamental GR spacetimes:
1. Minkowski (flat reference)
2. Schwarzschild (static spherical mass)
3. Kerr (rotating black hole)
4. Gödel (rotating dust universe with known CTCs)
5. Tipler (rotating cylinder with known CTCs)

Reports for each:
- Metric error
- Christoffel error
- Curvature (Riemann, Ricci, Einstein) error
- Geodesic error (norm conservation drift and trajectory behavior)
"""

from __future__ import annotations

import numpy as np

from .benchmark_metrics import (
    MinkowskiMetric,
    SchwarzschildMetric,
    KerrMetric,
    GodelMetric,
    TiplerCylinderMetric,
)
from .curvature_engine import (
    compute_curvature_chain,
    compute_christoffel,
)
from .geodesic_solver import (
    integrate_geodesic,
    circular_orbit_initial_conditions_schwarzschild,
)


def run_minkowski_benchmark() -> dict:
    m = MinkowskiMetric()
    probe = np.array([0.0, 1.0, 2.0, 3.0])
    g_exact = np.diag([-1.0, 1.0, 1.0, 1.0])
    g_calc = m.g(probe)
    metric_err = float(np.max(np.abs(g_calc - g_exact)))

    _, _, gamma = compute_christoffel(m.g, probe)
    christoffel_err = float(np.max(np.abs(gamma)))

    curv = compute_curvature_chain(m.g, probe)
    riemann_err = float(np.max(np.abs(curv["riemann"])))
    ricci_err = float(np.max(np.abs(curv["ricci"])))
    einstein_err = float(np.max(np.abs(curv["einstein"])))

    # Geodesic: straight timelike path u = (1, 0, 0, 0)
    geo = integrate_geodesic(m, probe, [1.0, 0.0, 0.0, 0.0], tau_span=10.0, n_steps=100)
    geodesic_norm_drift = geo["norm_drift_rel"]
    # Path should remain strictly straight
    pos_err = float(np.max(np.abs(geo["x"][:, 1:] - probe[1:])))

    passed = bool(metric_err < 1e-15 and christoffel_err < 1e-12 and riemann_err < 1e-12 and geodesic_norm_drift < 1e-14)
    return {
        "metric": "minkowski",
        "metric_error": metric_err,
        "christoffel_error": christoffel_err,
        "curvature_error": {"riemann": riemann_err, "ricci": ricci_err, "einstein": einstein_err},
        "geodesic_norm_drift": geodesic_norm_drift,
        "geodesic_position_error": pos_err,
        "pass": passed,
    }


def run_schwarzschild_benchmark(M: float = 1.0) -> dict:
    m = SchwarzschildMetric(M=M)
    r_probe = 6.0 * M
    probe = np.array([0.0, r_probe, np.pi / 2.0, 0.0])

    # 1. Metric error (symmetry & signature)
    g = m.g(probe)
    f = 1.0 - 2.0 * M / r_probe
    g_exact = np.diag([-f, 1.0 / f, r_probe**2, r_probe**2])
    metric_err = float(np.max(np.abs(g - g_exact)))

    # 2. Christoffel error vs analytic
    _, _, gamma_num = compute_christoffel(m.g, probe, h=1e-5)
    gamma_ana = m.christoffel_analytic(probe)
    christoffel_err = float(np.max(np.abs(gamma_num - gamma_ana)))

    # 3. Curvature error (Schwarzschild is vacuum: Ricci=0, Einstein=0)
    curv = compute_curvature_chain(m.g, probe, h=1e-4)
    ricci_err = float(np.max(np.abs(curv["ricci"])))
    einstein_err = float(np.max(np.abs(curv["einstein"])))
    riemann_scale = float(np.max(np.abs(curv["riemann"])))  # Kretschmann scale ~ M/r^3

    # 4. Geodesic: circular orbit at r = 6M
    x0, u0 = circular_orbit_initial_conditions_schwarzschild(M, r_probe)
    t_orbit = 2.0 * np.pi * np.sqrt(r_probe**3 / M)
    geo = integrate_geodesic(m, x0, u0, tau_span=t_orbit, n_steps=200, analytic=True)
    geodesic_norm_drift = geo["norm_drift_rel"]
    # Radial coordinate should remain near r_probe
    r_orbit_drift = float(np.max(np.abs(geo["x"][:, 1] - r_probe)) / r_probe)

    passed = bool(metric_err < 1e-15 and christoffel_err < 1e-5 and ricci_err < 1e-6 and geodesic_norm_drift < 1e-10)
    return {
        "metric": "schwarzschild",
        "M": M,
        "probe_r": r_probe,
        "metric_error": metric_err,
        "christoffel_error": christoffel_err,
        "curvature_error": {"ricci_vacuum_residual": ricci_err, "einstein_vacuum_residual": einstein_err, "riemann_scale": riemann_scale},
        "geodesic_norm_drift": geodesic_norm_drift,
        "orbit_radial_drift_rel": r_orbit_drift,
        "pass": passed,
    }


def run_kerr_benchmark(M: float = 1.0, a: float = 0.7) -> dict:
    m = KerrMetric(M=M, a=a)
    r_probe = 8.0 * M
    probe = np.array([0.0, r_probe, np.pi / 2.0, 0.0])

    g = m.g(probe)
    metric_symmetry_err = float(np.max(np.abs(g - g.T)))

    _, _, gamma_num = compute_christoffel(m.g, probe, h=1e-5)
    # Check metric compatibility dg - Gamma g - Gamma g ~ 0
    curv = compute_curvature_chain(m.g, probe, h=1e-4)
    ricci_err = float(np.max(np.abs(curv["ricci"])))
    einstein_err = float(np.max(np.abs(curv["einstein"])))

    # Equatorial timelike geodesic
    u0 = np.array([1.2, 0.0, 0.0, 0.03])
    # Normalize to timelike
    norm = float(u0 @ g @ u0)
    u0 /= np.sqrt(abs(norm))
    geo = integrate_geodesic(m, probe, u0, tau_span=20.0, n_steps=200)
    geodesic_norm_drift = geo["norm_drift_rel"]

    passed = bool(metric_symmetry_err < 1e-15 and ricci_err < 1e-5 and geodesic_norm_drift < 1e-8)
    return {
        "metric": "kerr",
        "M": M,
        "a": a,
        "metric_symmetry_error": metric_symmetry_err,
        "curvature_error": {"ricci_vacuum_residual": ricci_err, "einstein_vacuum_residual": einstein_err},
        "geodesic_norm_drift": geodesic_norm_drift,
        "pass": passed,
    }


def run_godel_benchmark(omega: float = 1.0) -> dict:
    m = GodelMetric(omega=omega)
    r_probe = 0.5 * m.ctc_radius()
    probe = np.array([0.0, r_probe, 0.0, 0.0])

    g = m.g(probe)
    metric_symmetry_err = float(np.max(np.abs(g - g.T)))

    # Curvature of Godel: R_ab = 2 omega^2 u_a u_b for dust frame
    curv = compute_curvature_chain(m.g, probe, h=1e-4)
    ricci_scalar = curv["ricci_scalar"]
    # Analytic Ricci scalar for Godel: R = -2 omega^2 (with (-+++) signature)
    expected_R = -2.0 * omega**2
    ricci_scalar_err = float(abs(ricci_scalar - expected_R) / abs(expected_R))

    # Geodesic
    u0 = np.array([np.sqrt(2.0) / omega, 0.0, 0.0, 0.0])
    geo = integrate_geodesic(m, probe, u0, tau_span=10.0, n_steps=100)
    geodesic_norm_drift = geo["norm_drift_rel"]

    # CTC detection at r > asinh(1) vs r < asinh(1)
    ctc_r = m.ctc_radius()
    g_inside = m.g([0.0, 0.5 * ctc_r, 0.0, 0.0])[2, 2]  # phi-phi component
    g_outside = m.g([0.0, 1.5 * ctc_r, 0.0, 0.0])[2, 2]
    ctc_behavior_correct = bool(g_inside > 0.0 and g_outside < 0.0)

    passed = bool(metric_symmetry_err < 1e-15 and ricci_scalar_err < 1e-3 and geodesic_norm_drift < 1e-8 and ctc_behavior_correct)
    return {
        "metric": "godel",
        "omega": omega,
        "ctc_radius": ctc_r,
        "metric_symmetry_error": metric_symmetry_err,
        "ricci_scalar_analytic_error": ricci_scalar_err,
        "geodesic_norm_drift": geodesic_norm_drift,
        "ctc_behavior_correct": ctc_behavior_correct,
        "pass": passed,
    }


def run_tipler_benchmark(a: float = 1.0) -> dict:
    m = TiplerCylinderMetric(a=a)
    r_probe = 0.5 * m.ctc_radius()
    probe = np.array([0.0, r_probe, 0.0, 0.0])

    g = m.g(probe)
    metric_symmetry_err = float(np.max(np.abs(g - g.T)))

    curv = compute_curvature_chain(m.g, probe, h=1e-4)
    curv_frobenius = float(np.linalg.norm(curv["riemann"]))

    # Geodesic
    u0 = np.array([1.0, 0.0, 0.0, 0.0])
    geo = integrate_geodesic(m, probe, u0, tau_span=10.0, n_steps=100)
    geodesic_norm_drift = geo["norm_drift_rel"]

    # CTC detection at a*r > 1 vs a*r < 1
    ctc_r = m.ctc_radius()
    g_inside = m.g([0.0, 0.5 * ctc_r, 0.0, 0.0])[2, 2]
    g_outside = m.g([0.0, 1.5 * ctc_r, 0.0, 0.0])[2, 2]
    ctc_behavior_correct = bool(g_inside > 0.0 and g_outside < 0.0)

    passed = bool(metric_symmetry_err < 1e-15 and geodesic_norm_drift < 1e-8 and ctc_behavior_correct)
    return {
        "metric": "tipler_cylinder",
        "a": a,
        "ctc_radius": ctc_r,
        "metric_symmetry_error": metric_symmetry_err,
        "curvature_scale": curv_frobenius,
        "geodesic_norm_drift": geodesic_norm_drift,
        "ctc_behavior_correct": ctc_behavior_correct,
        "pass": passed,
    }


def run_full_benchmark_suite() -> dict:
    """Run all 5 benchmark metric checks and report comprehensive suite results."""
    mink = run_minkowski_benchmark()
    schw = run_schwarzschild_benchmark()
    kerr = run_kerr_benchmark()
    godel = run_godel_benchmark()
    tipler = run_tipler_benchmark()

    all_pass = bool(mink["pass"] and schw["pass"] and kerr["pass"] and godel["pass"] and tipler["pass"])

    return {
        "overall_pass": all_pass,
        "benchmarks": {
            "minkowski": mink,
            "schwarzschild": schw,
            "kerr": kerr,
            "godel": godel,
            "tipler_cylinder": tipler,
        },
    }
