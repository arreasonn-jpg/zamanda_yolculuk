"""Causal-structure analyzer (ZT-006.3, PHASE 5/6 groundwork).

Classifies intervals and curves with the project-wide convention

    ds^2 < 0  -> timelike
    ds^2 = 0  -> null
    ds^2 > 0  -> spacelike

(signature (-,+,+,+), consistent with docs/physics.md).

No CTC conclusion is produced by this module alone; it only classifies
what is given to it.
"""

from __future__ import annotations

import numpy as np

TIMELIKE = "timelike"
NULL = "null"
SPACELIKE = "spacelike"


def classify_interval(ds2, tol=1e-12):
    """Classify a single line element value."""
    if ds2 < -tol:
        return TIMELIKE
    if ds2 > tol:
        return SPACELIKE
    return NULL


def classify_curve(metric, points, periodic_coord=None, period=2 * np.pi):
    """Classify a piecewise curve given by coordinate points (N, 4).

    The segment between consecutive points is evaluated with the metric at
    the segment midpoint.  If ``periodic_coord`` is given, coordinate
    differences along that index are wrapped into [-period/2, period/2]
    before computing dx (needed for phi-type angular coordinates).

    Returns a dict with per-segment ds^2, per-segment class, aggregate
    classification, and closure diagnostics.
    """
    pts = np.asarray(points, float)
    if pts.ndim != 2 or pts.shape[1] != 4:
        raise ValueError("points must have shape (N, 4)")
    n = len(pts)
    if n < 2:
        raise ValueError("need at least 2 points")

    ds2 = np.zeros(n - 1)
    for i in range(n - 1):
        dx = pts[i + 1] - pts[i]
        if periodic_coord is not None:
            dx[periodic_coord] = _wrap(dx[periodic_coord], period)
        mid = 0.5 * (pts[i + 1] + pts[i])
        if periodic_coord is not None:
            mid[periodic_coord] = _wrap(mid[periodic_coord], period)
        ds2[i] = metric.ds2(mid, dx)

    classes = [classify_interval(v) for v in ds2]
    if all(c == TIMELIKE for c in classes):
        overall = TIMELIKE
    elif all(c == SPACELIKE for c in classes):
        overall = SPACELIKE
    elif all(c == NULL for c in classes):
        overall = NULL
    elif all(c in (TIMELIKE, NULL) for c in classes):
        overall = "timelike_or_null"
    else:
        overall = "mixed"

    # closure in the covering space (angular coordinate unwrapped)
    d_close = pts[-1].copy()
    d_close -= pts[0]
    if periodic_coord is not None:
        d_close[periodic_coord] = _wrap(d_close[periodic_coord], period)

    return {
        "ds2": ds2,
        "segment_classes": classes,
        "overall": overall,
        "closed_covering_space": bool(np.allclose(d_close, 0.0, atol=1e-9)),
        "proper_time": float(np.sum(np.sqrt(np.maximum(-ds2, 0.0)))),
    }


def _wrap(d, period):
    return (d + period / 2.0) % period - period / 2.0
