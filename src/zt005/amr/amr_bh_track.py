"""Track BH positions from Atilde magnitude (v12.12-alpha1)."""
from __future__ import annotations
import numpy as np


def find_bh_centers(s, L, N, axis=0):
    """Return (x_left, x_right) BH center positions along axis.

    Simplification: punctures are along x-axis. Sum |Atilde|^2 over y,z,
    1D profile along x. Find two largest peaks.
    """
    A2 = np.einsum('...ij,...ij->...', s.Atilde, s.Atilde)
    if axis == 0:
        profile = A2.sum(axis=(1, 2))
    elif axis == 1:
        profile = A2.sum(axis=(0, 2))
    else:
        profile = A2.sum(axis=(0, 1))
    x = np.linspace(-L/2, L/2, N)
    # Split into left/right halves; find argmax in each
    mid = N // 2
    left_peak = np.argmax(profile[:mid])
    right_peak = np.argmax(profile[mid:]) + mid
    return x[left_peak], x[right_peak], profile


def pair_separation(x_left, x_right):
    return abs(x_right - x_left)


def peak_centroids(s, L, frac=0.05, radius=2):
    """2D BH centers (x, y) for an orbital pair, in grid coordinates.

    v14.0: an orbiting pair leaves the pure x-axis. Meter rules
    (lessons from the v14.0 step-0 calibration):
    - strict half-space split (x < 0 / x > 0): the exact-center node
      of an odd grid belongs to NEITHER side, or it tugs one side's
      centroid inward (an asymmetric bias with symmetric data);
    - per-side acceptance is limited to a local box of `radius` cells
      around the side's |Atilde|^2 argmax, so the shared inter-puncture
      field cannot tilt the weighted centroid;
    - the centroid is |Atilde|^2-weighted inside that box.

    Returns (c_left_xy, c_right_xy, A2_field).
    """
    A2 = np.einsum('...ij,...ij->...', s.Atilde, s.Atilde)
    N = A2.shape[0]
    x = np.linspace(-L / 2, L / 2, N)
    all_cells = np.argwhere(np.ones(A2.shape, dtype=bool))
    cx = x[all_cells[:, 0]]
    out = []
    for left in (True, False):
        side = (cx < 0.0) if left else (cx > 0.0)
        side_cells = all_cells[side]
        vals = A2[side_cells[:, 0], side_cells[:, 1], side_cells[:, 2]]
        pk = side_cells[int(np.argmax(vals))]
        box = (np.abs(side_cells[:, 0] - pk[0]) <= radius) & \
              (np.abs(side_cells[:, 1] - pk[1]) <= radius) & \
              (np.abs(side_cells[:, 2] - pk[2]) <= radius)
        c = side_cells[box]
        w = A2[c[:, 0], c[:, 1], c[:, 2]]
        wsum = float(np.sum(w))
        px = float(np.sum(x[c[:, 0]] * w) / wsum)
        py = float(np.sum(x[c[:, 1]] * w) / wsum)
        out.append((px, py))
    c_left, c_right = out
    return c_left, c_right, A2


def alpha_pair(s, L):
    """Alpha-funnel meter: per-half low-lapse basin centroid + depth.

    v14.0 second, independent locator for the moving-puncture gauge:
    after the gauge transient the lapse collapses into funnels at the
    punctures. Calibration lesson (single-grid probe to t = 6): on the
    coarse grid (h = 1.25) no funnel forms within 6 M (alpha_min
    RISES 0.67 -> 0.87 and the basin wanders); the meter is expected
    to become meaningful only at finer resolution / later times. It
    is recorded alongside the Atilde meter so the long run can decide
    empirically which locator tracks the inspiral.

    Per half-space (x < 0 / x > 0): the basin is the set of cells
    with alpha below a cut 10% into the depth range; the centroid is
    weighted by 1 / (alpha - depth), biasing toward the deepest
    cells. Returns (c_left_xy, c_right_xy, min_depths).
    """
    al = s.alpha
    N = al.shape[0]
    x = np.linspace(-L / 2, L / 2, N)
    all_cells = np.argwhere(np.ones(al.shape, dtype=bool))
    cx = x[all_cells[:, 0]]
    out, depths = [], []
    for left in (True, False):
        side = (cx < 0.0) if left else (cx > 0.0)
        sc = all_cells[side]
        vals = al[sc[:, 0], sc[:, 1], sc[:, 2]]
        depth = float(vals.min())
        cut = depth + 0.1 * (vals.max() - depth)
        c = sc[vals < cut]
        w = 1.0 / (al[c[:, 0], c[:, 1], c[:, 2]] - depth + 1e-2)
        wsum = float(np.sum(w))
        px = float(np.sum(x[c[:, 0]] * w) / wsum)
        py = float(np.sum(x[c[:, 1]] * w) / wsum)
        out.append((px, py))
        depths.append(depth)
    c_left, c_right = out
    return c_left, c_right, tuple(depths)


__all__ = ["find_bh_centers", "pair_separation", "peak_centroids",
           "alpha_pair"]
