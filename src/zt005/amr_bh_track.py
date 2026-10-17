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


__all__ = ["find_bh_centers", "pair_separation"]
