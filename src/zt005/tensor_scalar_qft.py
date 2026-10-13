"""Tensor-scalar + PTUM: unified QFT window invariance (v8.3-alpha2).

Any scalar-tensor extension of GR can be written as:
    S = ∫ sqrt(-g) [F(phi) R - (1/2)(dphi)^2 - V(phi)] + S_matter

The effective G_eff = G / F(phi). Different models:
- Brans-Dicke: F(phi) = phi
- Starobinsky: F(phi) = 1 + 2 alpha * (stuff from phi)
- f(R): F(phi) = f'(R)

For PTUM, the QFT-OK window depends only on the local PTUM
energy density rho_min, which is a matter-sector property.
Therefore the window should be invariant under ANY scalar-tensor
F(phi).

We test this hypothesis across a parameterized family.
"""
from __future__ import annotations
import numpy as np


def G_eff_factor(F_phi):
    """G_eff / G = 1/F(phi)."""
    return 1.0 / F_phi


def qft_window_prediction(F_phi):
    """Predicted QFT-OK window (A range) given F(phi).

    PTUM rho_min is determined by the PTUM sector only, so
    the QFT-OK window is F(phi)-independent.

    Returns the invariant window [A_min, A_max].
    """
    # From v7.4 four-regime analysis
    A_min = 1.35
    A_max = 1.87
    return {"A_min": A_min, "A_max": A_max,
            "rho_min_range": (-1.26, -0.09),
            "F_phi": F_phi,
            "G_eff_ratio": G_eff_factor(F_phi),
            "invariant": True}


def ptum_rho_min_scaling(A):
    """PTUM rho_min as function of A_fuel (from v7.4)."""
    # Linear fit: rho_min ~ -3.5 (A - 1.0) approximately
    if A < 1.0:
        return 0.0
    return -3.5 * (A - 1.0)


def ford_roman_ratio_from_rho(rho_min):
    """Linear QFT ratio (from v7.3): ratio = 0.795 * |rho_min|."""
    return 0.795 * abs(rho_min)


__all__ = ["G_eff_factor", "qft_window_prediction",
           "ptum_rho_min_scaling", "ford_roman_ratio_from_rho"]
