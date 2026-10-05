"""Closed-timelike-curve (CTC) benchmark search (ZT-006.3, PHASE 6).

The first question any CTC claim must pass is:

    "Can the solver find a CTC in a spacetime where one is KNOWN to exist?"

This module implements that known-answer test.  It searches the simplest
analytic CTC family — coordinate circles (t = const, r = const, z = const,
phi: 0 -> 2*pi) — and reports whether the circle is closed and everywhere
timelike.  A curve is reported as a CTC candidate only when BOTH hold.

Known analytic answers used for validation:

* Gödel universe:        phi-circles timelike  <=>  r > asinh(1)
* Tipler cylinder (idealized): phi-circles timelike  <=>  a*r > 1
* Minkowski / Schwarzschild exterior / Kerr exterior / Morris-Thorne:
                         phi-circles are spacelike everywhere in the domain.

Important honesty constraints:

* Finding a CTC in a benchmark metric says NOTHING about the device model.
* Failing to find one in the device model (which is not even a full metric
  yet) says nothing either — the geodesic/CTC pipeline is validated only
  against these known answers at this stage.
"""

from __future__ import annotations

import numpy as np

from .causal_structure import TIMELIKE, classify_curve


def phi_circle_points(r, z_or_theta, n_points=129, t0=0.0,
                      coord_layout="trpz"):
    """Sample a constant-t, constant-r circle.

    coord_layout:
      "trpz"  -> coordinates (t, r, phi, z)   [Gödel, Tipler]
      "trtp"  -> coordinates (t, r, theta, phi) [spherical-style]
    """
    phi = np.linspace(0.0, 2.0 * np.pi, n_points)
    pts = np.zeros((n_points, 4))
    pts[:, 0] = t0
    pts[:, 1] = r
    if coord_layout == "trpz":
        pts[:, 2] = phi
        pts[:, 3] = z_or_theta
    elif coord_layout == "trtp":
        pts[:, 2] = z_or_theta   # theta
        pts[:, 3] = phi          # phi
    else:
        raise ValueError("unknown coord_layout")
    return pts


def check_phi_circle_ctc(metric, r, z_or_theta=0.0, n_points=129,
                         coord_layout="trpz"):
    """Evaluate the constant-t phi-circle as a CTC candidate."""
    pts = phi_circle_points(r, z_or_theta, n_points,
                            coord_layout=coord_layout)
    periodic = 2 if coord_layout == "trpz" else 3
    report = classify_curve(metric, pts, periodic_coord=periodic)
    report["radius"] = float(r)
    report["is_ctc_candidate"] = bool(
        report["closed_covering_space"] and report["overall"] == TIMELIKE)
    return report


def godel_ctc_benchmark(omega=1.0, n_points=129):
    """Known-answer test in Gödel spacetime.

    Returns dict with detection results below and above the analytic CTC
    radius r_G = asinh(1), plus the agreement flag.
    """
    from .benchmark_metrics import GodelMetric
    m = GodelMetric(omega=omega)
    r_in = 0.5 * m.ctc_radius()
    r_out = 1.5 * m.ctc_radius()
    inside = check_phi_circle_ctc(m, r_in, n_points=n_points,
                                  coord_layout="trpz")
    outside = check_phi_circle_ctc(m, r_out, n_points=n_points,
                                   coord_layout="trpz")
    return {
        "metric": "godel",
        "ctc_radius_analytic": m.ctc_radius(),
        "r_inside": r_in, "detected_inside": inside["is_ctc_candidate"],
        "r_outside": r_out, "detected_outside": outside["is_ctc_candidate"],
        "agrees_with_analytic": bool(
            (not inside["is_ctc_candidate"]) and outside["is_ctc_candidate"]),
        "inside": inside, "outside": outside,
    }


def tipler_ctc_benchmark(a=1.0, n_points=129):
    """Known-answer test in the idealized Tipler-cylinder metric."""
    from .benchmark_metrics import TiplerCylinderMetric
    m = TiplerCylinderMetric(a=a)
    r_in = 0.5 * m.ctc_radius()
    r_out = 1.5 * m.ctc_radius()
    inside = check_phi_circle_ctc(m, r_in, n_points=n_points,
                                  coord_layout="trpz")
    outside = check_phi_circle_ctc(m, r_out, n_points=n_points,
                                   coord_layout="trpz")
    return {
        "metric": "tipler_cylinder_idealized",
        "ctc_radius_analytic": m.ctc_radius(),
        "r_inside": r_in, "detected_inside": inside["is_ctc_candidate"],
        "r_outside": r_out, "detected_outside": outside["is_ctc_candidate"],
        "agrees_with_analytic": bool(
            (not inside["is_ctc_candidate"]) and outside["is_ctc_candidate"]),
        "inside": inside, "outside": outside,
    }


def negative_benchmarks(n_points=65):
    """Metrics where the phi-circle search must NOT find a CTC.

    Any detection here is a solver defect, not a physics result.
    """
    from .benchmark_metrics import (MinkowskiMetric, MorrisThorneMetric,
                                    SchwarzschildMetric)
    results = {}

    m = MinkowskiMetric()
    # circle of coordinate radius R in the (x, y) plane at t = 0
    phi = np.linspace(0, 2 * np.pi, n_points)
    pts = np.zeros((n_points, 4))
    pts[:, 1] = 1.0 * np.cos(phi)
    pts[:, 2] = 1.0 * np.sin(phi)
    rep = classify_curve(m, pts)
    results["minkowski"] = {"is_ctc_candidate": False,
                            "overall": rep["overall"],
                            "found": rep["overall"] == TIMELIKE}

    m = SchwarzschildMetric(M=1.0)
    rep = check_phi_circle_ctc(m, 6.0, z_or_theta=np.pi / 2,
                               n_points=n_points, coord_layout="trtp")
    results["schwarzschild_r6M"] = {
        "is_ctc_candidate": rep["is_ctc_candidate"],
        "overall": rep["overall"],
        "found": rep["is_ctc_candidate"]}

    m = MorrisThorneMetric(b0=1.0)
    # constant-l circle at the throat, theta = pi/2, coords (t,l,theta,phi)
    pts = np.zeros((n_points, 4))
    pts[:, 1] = 0.0
    pts[:, 2] = np.pi / 2
    pts[:, 3] = np.linspace(0, 2 * np.pi, n_points)
    rep = classify_curve(m, pts, periodic_coord=3)
    results["morris_thorne_throat"] = {
        "is_ctc_candidate": False,
        "overall": rep["overall"],
        "found": rep["overall"] == TIMELIKE}

    results["all_clean"] = bool(not any(r["found"] for r in results.values()
                                        if isinstance(r, dict)))
    return results
