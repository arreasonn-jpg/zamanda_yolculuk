"""CTC Search Engine & "1 Second Past" Backward Time-Travel Gate.

Addresses Section 12 and 13:
- CTC search:
    Finds closed loops x^mu(lambda) with x(0) == x(1) and ds^2 < 0 everywhere.
- Direct Criterion for "1 Second Past":
    Delta t = t_final - t_initial < 0  (specifically Delta t <= -1.0 s)
    WHILE proper time tau_traveler > 0 (traveler ages forward while moving backward in external coordinate time).
"""

from __future__ import annotations

import numpy as np

from .causal_structure import classify_curve, TIMELIKE
from .benchmark_metrics import GodelMetric, TiplerCylinderMetric, MinkowskiMetric

C = 299_792_458.0


def evaluate_closed_curve_ctc(
    metric,
    points: np.ndarray,
    tol_closure: float = 1e-4,
    periodic_coord: int | None = None,
) -> dict:
    """Evaluate whether a sampled parameterized curve constitutes a valid Closed Timelike Curve."""
    pts = np.asarray(points, float)
    if len(pts) < 3:
        raise ValueError("Need at least 3 points to evaluate a closed curve")

    # 1. Closure check
    d_close = pts[-1].copy() - pts[0]
    if periodic_coord is not None:
        d_close[periodic_coord] = (d_close[periodic_coord] + np.pi) % (2.0 * np.pi) - np.pi
    closure_dist = float(np.linalg.norm(d_close))
    is_closed = bool(closure_dist <= tol_closure)

    # 2. Causal classification along curve
    report = classify_curve(metric, pts, periodic_coord=periodic_coord)
    is_everywhere_timelike = bool(report["overall"] == TIMELIKE)

    # CTC Candidate flag: closed AND everywhere timelike
    is_ctc = bool(is_closed and is_everywhere_timelike)

    return {
        "is_ctc": is_ctc,
        "is_closed": is_closed,
        "closure_distance": closure_dist,
        "is_everywhere_timelike": is_everywhere_timelike,
        "overall_causal_class": report["overall"],
        "proper_time_tau": report["proper_time"],
        "n_segments": len(report["ds2"]),
    }


def check_one_second_backward_travel_criterion(
    metric,
    curve_points: np.ndarray,
    target_dt_seconds: float = -1.0,
    periodic_coord: int | None = None,
) -> dict:
    """Evaluate the explicit '1 Second Past' criterion.

    Requirements:
    1. Delta t = t_final - t_initial <= target_dt_seconds (< 0)
    2. Proper time tau_traveler > 0 (traveler's internal clock ticks forward)
    3. Entire worldline is physically timelike (ds^2 < 0)
    """
    pts = np.asarray(curve_points, float)
    t_initial = float(pts[0, 0] / C)
    t_final = float(pts[-1, 0] / C)
    delta_t = t_final - t_initial

    # Evaluate causal structure
    curve_eval = classify_curve(metric, pts, periodic_coord=periodic_coord)
    proper_time_tau = curve_eval["proper_time"]
    is_timelike = bool(curve_eval["overall"] == TIMELIKE)

    # Check the required criteria:
    reaches_past_time = bool(delta_t <= target_dt_seconds)
    proper_time_positive = bool(proper_time_tau > 0.0)
    criterion_met = bool(reaches_past_time and proper_time_positive and is_timelike)

    return {
        "one_second_backward_criterion_met": criterion_met,
        "target_dt_s": target_dt_seconds,
        "achieved_delta_t_s": delta_t,
        "proper_time_tau_s": proper_time_tau,
        "reaches_past_time": reaches_past_time,
        "proper_time_positive": proper_time_positive,
        "is_timelike_worldline": is_timelike,
        "initial_coordinate_time_s": t_initial,
        "final_coordinate_time_s": t_final,
    }


def scan_for_ctc_in_metric(
    metric,
    r_min: float,
    r_max: float,
    n_radii: int = 20,
    periodic_coord: int = 2,
) -> dict:
    """Scan radial domain of an axisymmetric metric for closed timelike curves."""
    radii = np.linspace(r_min, r_max, n_radii)
    ctc_candidates = []

    for r in radii:
        phi = np.linspace(0.0, 2.0 * np.pi, 129)
        # Circle at constant t=0, z=0
        pts = np.zeros((129, 4))
        pts[:, 0] = 0.0
        pts[:, 1] = r
        pts[:, 2] = phi
        pts[:, 3] = 0.0

        eval_res = evaluate_closed_curve_ctc(metric, pts, tol_closure=1e-3, periodic_coord=periodic_coord)
        if eval_res["is_ctc"]:
            ctc_candidates.append({
                "radius": float(r),
                "proper_time": eval_res["proper_time_tau"],
                "g_phiphi": float(metric.g(pts[0])[2, 2]),
            })

    has_ctc = len(ctc_candidates) > 0
    min_ctc_radius = min((c["radius"] for c in ctc_candidates), default=None)

    return {
        "metric_name": getattr(metric, "name", "custom"),
        "has_ctc": has_ctc,
        "min_ctc_radius": min_ctc_radius,
        "n_ctc_detected": len(ctc_candidates),
        "candidates": ctc_candidates,
    }
